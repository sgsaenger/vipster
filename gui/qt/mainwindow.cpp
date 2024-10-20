#include "version.h"
#include "mainwindow.h"
#include "ui_viewport.h"
#include "savefmtdialog.h"
#include "mainwidgets.h"
#include "toolwidgets.h"
#include "glwidget.h"

#include <QDockWidget>
#include <QMessageBox>
#include <QMenuBar>
#include <QFileDialog>
#include <QInputDialog>
#include <QApplication>
#include <QStyle>
#include <QToolBar>

#include <filesystem>
#include <fmt/format.h>

using namespace Vipster;
namespace fs = std::filesystem;

MainWindow * MainWindow::_instance = nullptr;

MainWindow::MainWindow(Vipster::ConfigState state, QDir path):
    QMainWindow{},
    _config{std::move(state)},
    path{path}
{
    if (_instance) {
        throw Error{"Only one instance of MainWindow must exist"};
    } else {
        _instance = this;
    }

    // set up GUI
    setWindowIcon(QIcon{":/images/vipster.png"});
    setupFileMenu();
    setupEditMenu();
    setupHelpMenu();
    setupMainWidgets();
    setupToolWidgets();
    setupViewports();

    // global event handlers

    // Ensure bonds and overlap information is up to date
    auto updateBonds = [&](const Step &step){
        if (getState(step).automatic_bonds ||
            _config.settings.overlap.val) {
            const_cast<Step&>(step).generateBonds(!getState(step).automatic_bonds);
        }
    };
    connect(this, &MainWindow::stepChanged, this, updateBonds);
}

/********************************
 * Operations on molecule files *
 ********************************/
void MainWindow::readFile()
{
    // File dialog
    QFileDialog fileDiag{this};
    fileDiag.setDirectory(path);
    // TODO: limited to one file for now
    fileDiag.setFileMode(QFileDialog::ExistingFile);
    // Format dialog
    QStringList formats{};
    for(const auto plug: _config.plugins){
        if(plug->parser)
            formats << QString::fromStdString(plug->name);
    }
    if(fileDiag.exec() != 0){
        auto files = fileDiag.selectedFiles();
        path = fileDiag.directory();
        // if cancelled by user, just return
        if(files.empty()) return;
        // else try to read file
        const auto& file = files[0].toStdString();
        // guess format or request from user
        auto plugin = guessFmt(file, _config.plugins);
        if (!plugin){
            bool got_fmt{false};
            auto fmt_s = QInputDialog::getItem(this, "Select format", "Format:",
                                               formats, 0, false, &got_fmt);
            // if the user selected the format, read the file
            if(got_fmt){
                auto fmt = std::find_if(_config.plugins.begin(), _config.plugins.end(),
                    [&](const auto& plug){return plug->name.c_str() == fmt_s;});
                if(fmt == _config.plugins.end())
                    throw Error{"Invalid format in loadMol occured"};
                plugin = *fmt;
            }
        }
        try {
            if(plugin){
                // try to open file
                newIOData(Vipster::readFile(file, plugin));
            }
        }catch (const IOError &e){
            QMessageBox msg{this};
            msg.setText(QString{"Could not open file \""}+file.c_str()+"\":\n"+e.what());
            msg.exec();
        }catch (const Error &e){
            QMessageBox msg{this};
            msg.setText(QString{"Could not open file \""}+file.c_str()+"\":\n"+e.what());
            msg.exec();
        }
    }
}

void MainWindow::saveFile()
{
    QFileDialog fileDiag{this};
    fileDiag.setDirectory(path);
    fileDiag.setAcceptMode(QFileDialog::AcceptSave);
    if(fileDiag.exec() == QDialog::Accepted){
        auto target = fileDiag.selectedFiles()[0].toStdString();
        path = fileDiag.directory();
        SaveFmtDialog sfd{_config.plugins, this};
        if(sfd.exec() == QDialog::Accepted){
            try{
                writeFile(target, sfd.plugin, curMol(),
                          curVP->moldata[&curMol()].curStep-1,
                          sfd.getParam(), sfd.getPreset());
            }catch(const IOError& e){
                QMessageBox msg{this};
                msg.setText(QString{"Could not write file \""}+target.c_str()+"\":\n"+e.what());
                msg.exec();
            }
        }
    }
}

void MainWindow::newIOData(IOTuple &&t)
{
    newMol(std::move(std::get<Molecule>(t)));

    const auto& name = molecules.back().name;
    if (auto &param = std::get<std::optional<Parameter>>(t)) {
        paramWidget->registerParam(name, std::move(*param));
        emit parameterListChanged();
    }

    for (auto &extra: std::get<DataList>(t)) {
        data.emplace_back(std::move(extra));
        emit dataListChanged();
    }
}


/****************************
 * Molecule & Step handling *
 ****************************/
void MainWindow::newMol(Molecule &&mol){
    molecules.emplace_back(std::move(mol));
    molecules.back().getPTE().root = &_config.periodicTable;

    emit molListChanged(molecules);
}

Molecule& MainWindow::curMol()
{
    return *pCurMol;
}

void MainWindow::setActiveMol(Molecule &m)
{
    pCurMol = &m;
    emit activeMolChanged(m);
}

const Step& MainWindow::curStep()
{
    return *pCurStep;
}

const Step::selection& MainWindow::curSel()
{
    return *pCurSel;
}

void MainWindow::setActiveStep(Step &step, Step::selection &sel)
{
    pCurStep = &step;
    pCurSel = &sel;

    getState(step);
    emit activeStepChanged(step, sel);
}

void MainWindow::updateSelection(SelectionFilter filter)
{
    *pCurSel = pCurStep->select(filter);
    emit selChanged(*pCurSel);
}

MainWindow::StepState& MainWindow::getState(const Step &step)
{
    // initialize state if required
    if (stepdata.find(&step) == stepdata.end()) {
        stepdata.emplace(&step, StepState{
            // when bonds are present (e.g. in file), use manual mode
            step.getBonds().empty(),
            // default formatter is passthrough
            const_cast<Step&>(step).asFmt(step.getFmt())
        });
    }
    return stepdata.at(&step);
}

void MainWindow::selectionToCopy(){
    copyBuf = *pCurSel;
    emit copyBufChanged(copyBuf);
}


/* Main widgets
 * Core functionality (Molecule-Properties, Auxiliary data)
 * Always open, stacked in left dock-area
 */
void MainWindow::setupMainWidgets()
{
    setDockOptions(QMainWindow::AllowTabbedDocks|QMainWindow::AnimatedDocks|QMainWindow::VerticalTabs);
#ifdef Q_OS_MACOS
    setDockOptions(dockOptions()^VerticalTabs);
#endif

    QDockWidget* firstDock{nullptr};
    for(const auto& pair: makeMainWidgets(this)){
        auto *tmp = new QDockWidget(this);
        tmp->setWidget(pair.first);
        tmp->setAllowedAreas(Qt::LeftDockWidgetArea);
        tmp->setFeatures(QDockWidget::DockWidgetMovable |
                         QDockWidget::DockWidgetFloatable);
        tmp->setWindowTitle(pair.second);
        addDockWidget(Qt::LeftDockWidgetArea, tmp);
        if(!firstDock){
            firstDock = tmp;
        }else{
            tabifyDockWidget(firstDock, tmp);
        }
        auto p = dynamic_cast<ParamWidget*>(pair.first);
        if (p) paramWidget = p;
        auto c = dynamic_cast<PresetWidget*>(pair.first);
        if (c) presetWidget = c;
    }
    firstDock->raise();
}

/* Tool widgets
 * Optional widgets with auxiliary functionality
 * Can be docked everywhere but left, opened&closed
 * Accessed via toolbar
 */
void MainWindow::setupToolWidgets()
{
    auto &toolBar = *addToolBar("toolBar");
    for(const auto& pair: makeToolWidgets(this)){
        auto *tmp = new QDockWidget(this);
        tmp->setWidget(pair.first);
        tmp->setAllowedAreas(Qt::BottomDockWidgetArea|
                             Qt::RightDockWidgetArea|
                             Qt::TopDockWidgetArea);
        tmp->setWindowTitle(pair.second);
        addDockWidget(Qt::RightDockWidgetArea, tmp);
        auto action = toolBar.addAction(pair.second);
        action->setCheckable(true);
        connect(action, &QAction::toggled, tmp, &QWidget::setVisible);
        tmp->hide();
    }
}

/* File menu
 * Load/save molecular and auxiliary data
 */
void MainWindow::setupFileMenu()
{
    auto &fileMenu = *menuBar()->addMenu("&File");

    // Create an empty, new molecule
    fileMenu.addAction("&New molecule",
                       [this](){newMol({});},
                       QKeySequence::New);

    // create a copy of existing data
    auto &copyMenu = *fileMenu.addMenu("&From existing molecule");
    copyMenu.addAction("Copy single &step",
                       [this](){
                           auto tmpMol = Molecule{curStep()};
                           tmpMol.name += " (copy of step " +
                                          std::to_string(curVP->moldata[&curMol()].curStep) + ")";
                           newMol(std::move(tmpMol));
                       });
    copyMenu.addAction("Copy &current selection",
                       [this](){
                           auto tmpMol = Molecule{curSel()};
                           tmpMol.name += " (copy of selection of step " +
                                          std::to_string(curVP->moldata[&curMol()].curStep) + ")";
                           newMol(std::move(tmpMol));
                       });
    copyMenu.addAction("Copy &trajectory",
                       [this](){
                           Molecule tmpMol = curMol();
                           tmpMol.name += " (copy)";
                           newMol(std::move(tmpMol));
                       });

    // load molecular data from file
    fileMenu.addAction("&Load molecule",
                       &MainWindow::readFile,
                       QKeySequence::Open);
    // store molecular data to file
    fileMenu.addAction("&Save molecule",
                       &MainWindow::saveFile,
                       QKeySequence::Save);

    // Separator
    fileMenu.addSeparator();
    // Create a nested menu that exposes all parameter presets per filetype plugin
    auto &paramMenu = *fileMenu.addMenu("Load parameter set");
    auto &presetMenu = *fileMenu.addMenu("Load IO preset");
    for(const auto& plug: config().plugins){
        if(plug->makeParam){
            // create a sub-menu for each param-enabled plugin
            auto* plug_menu = paramMenu.addMenu(QString::fromStdString(plug->name));
            // register each parameter preset as a separate action
            const auto& param_map = _config.parameters.at(plug);
            for(const auto& p: param_map){
                plug_menu->addAction(QString::fromStdString(p.first),
                                     this, &MainWindow::loadParam);
            }
        }
        if(plug->makePreset){
            // create a sub-menu for each IO-preset enabled plugin
            auto* plug_menu = presetMenu.addMenu(QString::fromStdString(plug->name));
            // register each IO preset as a separate action
            const auto& conf_map = _config.presets.at(plug);
            for(const auto& p: conf_map){
                plug_menu->addAction(QString::fromStdString(p.first),
                                     this, &MainWindow::loadPreset);
            }
        }
    }

    // Separator
    fileMenu.addSeparator();
    fileMenu.addAction("Screenshot (current Step)",
                       [this](){saveScreenshot();},
                       QKeySequence::Print);
    fileMenu.addAction("Screenshot (trajectory)",
                       [this](){saveScreenshots();});

    // Separator
    fileMenu.addSeparator();
    fileMenu.addAction("Exit Vipster",
                       [this](){close();},
                       Qt::ControlModifier|Qt::Key_Q);
}

/* Edit menu
 * Modify the currently selected Molecule/Step/Atoms
 */
void MainWindow::setupEditMenu()
{
    auto &editMenu = *menuBar()->addMenu("&Edit");

    // Create a new Atom
    auto *newAction = editMenu.addAction("&New atom",
        [this](){
            invokeOnStep(static_cast<void(Step::*)(const std::string &, const Vec&, const AtomProperties&)>(&Step::newAtom),
                          "C", Vec{}, AtomProperties{});
        },
        Qt::Key_N);

    // Delete selected Atom(s)
    auto delAtoms = [](Step &s, const Step::const_selection &sel)
    {
        std::set<size_t> indices{};
        for(const auto [idx, _]: sel.getAtoms().indices){
            indices.insert(idx);
        }
        for(auto it = indices.rbegin(); it != indices.rend(); ++it)
        {
            if(*it < s.getNat()){
                s.delAtom(*it);
            }
        }
    };
    auto *delAction = editMenu.addAction("&Delete atom(s)",
        [&](){
            invokeOnStep(delAtoms, curSel());
            updateSelection({});
        },
        Qt::Key_Delete);
    delAction->setEnabled(false);
    connect(this, &MainWindow::selChanged,
            delAction, [&](const Step::selection &sel){delAction->setEnabled(sel.getNat() > 0);});

    // Cut atoms
    auto *cutAction = editMenu.addAction("C&ut atom(s)",
        [&](){
            selectionToCopy();
            invokeOnStep(delAtoms, curSel());
            updateSelection({});
        },
        QKeySequence::Cut);
    cutAction->setEnabled(false);
    connect(this, &MainWindow::selChanged,
            cutAction, [&](const Step::selection &sel){cutAction->setEnabled(sel.getNat() > 0);});

    // Copy atoms
    auto *copyAction = editMenu.addAction("&Copy atom(s)",
        [&](){
            selectionToCopy();
        },
        QKeySequence::Copy);
    copyAction->setEnabled(false);
    connect(this, &MainWindow::selChanged,
            copyAction, [&](const Step::selection &sel){copyAction->setEnabled(sel.getNat() > 0);});

    // Paste atoms
    auto *pasteAction = editMenu.addAction("&Paste atom(s)",
        [&](){
            invokeOnStep([&](Step &s){
                if (copyBuf.getNat() > 0) {
                    s.newAtoms(copyBuf);
                }
            });
        },
        QKeySequence::Paste);
    pasteAction->setEnabled(false);
    connect(this, &MainWindow::copyBufChanged,
            pasteAction, [&](const std::optional<Step> &buf){pasteAction->setEnabled(buf->getNat() > 0);});

    // Separator
    editMenu.addSeparator();

    // Rename
    auto *renameAction = editMenu.addAction("&Rename atom(s)",
        [this](){
            auto newName = QInputDialog::getText(this, "Rename atoms",
                                             "Enter new Atom-type for selected atoms:")
                       .toStdString();
            invokeOnSelection([](Step::selection &sel, const std::string &name){
                for(auto& at: sel){
                    at.name = name;
                }
            }, newName);
        });
    connect(this, &MainWindow::selChanged,
            renameAction, [=](const Step::selection &sel){renameAction->setEnabled(sel.getNat() > 0);});

    // Hide
    auto *hideAction = editMenu.addAction("&Hide atom(s)",
          [&](){
              invokeOnSelection([](Step::selection &sel){
                  for(auto& at: sel){
                      at.properties->flags[AtomProperties::Hidden] = true;
                  }
              });
          });
      connect(this, &MainWindow::selChanged,
              hideAction, [&](const Step::selection &sel){hideAction->setEnabled(sel.getNat() > 0);});

    // Show
    auto *showAction = editMenu.addAction("&Show atom(s)",
        [&](){
            invokeOnSelection([](Step::selection &sel){
                for(auto& at: sel){
                    at.properties->flags[AtomProperties::Hidden] = false;
                }
            });
        });
    connect(this, &MainWindow::selChanged,
            showAction, [&](const Step::selection &sel){showAction->setEnabled(sel.getNat() > 0);});
}

void MainWindow::setupHelpMenu()
{
    auto &helpMenu = *menuBar()->addMenu("&Help");
    helpMenu.addAction("About Vipster", [&](){
        QMessageBox::about(this,
            QString("About Vipster"),
            QString("<h1>Vipster v" VIPSTER_VERSION "</h1>"
                    "<h3>(" VIPSTER_PLATFORM " " VIPSTER_ARCH ")</h3>"
                    "<p>"
                    "©Sebastian Gsänger, 2014 - " VIPSTER_YEAR
                    "<br>"
                    "<a href='https://sgsaenger.github.io/vipster'>Homepage</a>"
                    "<br>"
                    "<a href='https://github.com/sgsaenger/vipster'>Source</a>"
                    "</p>"
                    "<p>"
                    "This program is provided under the GPLv3."
                    "<br>"
                    "It uses<br>"
                    "<a href='https://github.com/nlohmann/json'>JSON for Modern C++</a>,<br>"
                    "<a href='https://github.com/CLIUtils/CLI11'>CLI11</a>,<br>"
                    "<a href='https://github.com/codeplea/tinyexpr'>TinyExpr</a>,<br>"
                    "<a href='https://github.com/fmtlib/fmt'>{fmt}</a>,<br>"
                    "<a href='https://github.com/catchorg/catch2'>Catch2</a><br>"
                    "and <a href='https://github.com/pybind/pybind11'>pybind11</a>."
                    "</p>"));
    });
    helpMenu.addAction("About Qt", &QApplication::aboutQt);
}

void MainWindow::setupViewports()
{
    // create first level of splitters and main viewport
    viewports.push_back(new ViewPort{this});
    vsplit = new QSplitter{this};
    vsplit->setOrientation(Qt::Vertical);
    vsplit->setChildrenCollapsible(false);
    hsplits.push_back(new QSplitter{vsplit});
    hsplits.back()->setChildrenCollapsible(false);
    hsplits.back()->addWidget(viewports.front());
    setCentralWidget(vsplit);
    curVP = viewports.front();
    curVP->setFrameShadow(QFrame::Shadow::Sunken);
    curVP->ui->closeButton->setDisabled(true);
}

void MainWindow::setActiveViewport(ViewPort* sender) {
    // deactivate other viewports
    for(auto& vp: viewports){
        vp->setFrameShadow(QFrame::Shadow::Raised);
    }
    sender->setFrameShadow(QFrame::Shadow::Sunken);
    // make sender current viewport
    curVP = sender;
    setActiveMol(*curVP->curMol);
    setActiveStep(*curVP->curStep, *curVP->curSel);
}

void MainWindow::splitViewportHoriz(ViewPort* sender) {
    viewports.front()->ui->closeButton->setEnabled(true);
    viewports.push_back(new ViewPort{*sender});
    for(auto& h: hsplits){
        for(auto& c: h->children()){
            if (c == sender){
                h->addWidget(viewports.back());
                return;
            }
        }
    }
    throw Error{"Could not determine horizontal splitter for viewport"};
}

void MainWindow::splitViewportVert(ViewPort* sender) {
    viewports.front()->ui->closeButton->setEnabled(true);
    viewports.push_back(new ViewPort{*sender});
    hsplits.push_back(new QSplitter{vsplit});
    hsplits.back()->setChildrenCollapsible(false);
    hsplits.back()->addWidget(viewports.back());
    vsplit->addWidget(hsplits.back());
}

void MainWindow::closeViewport(ViewPort* sender) {
    bool wasActive = sender == curVP;
    auto pos = std::find(viewports.begin(), viewports.end(), sender);
    if(pos == viewports.end()){
        throw Error{"Invalid viewport"};
    }else{
        viewports.erase(pos);
        delete sender;
        // clean-up hsplits if sender was last child
        for(auto it = hsplits.begin(); it != hsplits.end(); ++it){
            if(!(*it)->count()){
                delete *it;
                hsplits.erase(it);
                break;
            }
        }
    }
    // make sure we don't close last viewport
    if(viewports.size() == 1)
        viewports.front()->ui->closeButton->setDisabled(true);
    if(wasActive)
        viewports.front()->openGLWidget->setFocus();
}

const decltype (ParamWidget::params)& MainWindow::getParams() const noexcept
{
    return paramWidget->params;
}

const decltype (PresetWidget::presets)& MainWindow::getPresets() const noexcept
{
    return presetWidget->presets;
}

void MainWindow::loadParam()
{
    auto* s = static_cast<QAction*>(sender());
    auto* p = static_cast<QMenu*>(s->parent());
    auto fmt = std::find_if(_config.plugins.begin(), _config.plugins.end(),
                 [&](const auto& plug){return plug->name.c_str() == p->title();});
    if(fmt == _config.plugins.end()){
        throw Error{"Invalid parameter set"};
    }
    auto name = s->text().toStdString();
    auto pos = _config.parameters.at(*fmt).find(name);
    if(pos != _config.parameters.at(*fmt).end()){
        paramWidget->registerParam(name, pos->second);
    }else{
        throw Error("Invalid parameter set");
    }
}

void MainWindow::loadPreset()
{
    auto* s = static_cast<QAction*>(sender());
    auto* p = static_cast<QMenu*>(s->parent());
    auto fmt = std::find_if(_config.plugins.begin(), _config.plugins.end(),
                 [&](const auto& plug){return plug->name.c_str() == p->title();});
    if(fmt == _config.plugins.end()){
        throw Error{"Invalid IO-preset"};
    }
    auto name = s->text().toStdString();
    auto pos = _config.presets.at(*fmt).find(name);
    if(pos != _config.presets.at(*fmt).end()){
        presetWidget->registerPreset(name, pos->second);
    }else{
        throw Error("Invalid IO preset");
    }
}

void MainWindow::saveParam()
{
    if(!paramWidget->curParam){
        return;
    }
    bool ok;
    const auto &curParam = paramWidget->curParam;
    auto name = QInputDialog::getText(this, "Save parameter set", "Name of preset",
                                      QLineEdit::Normal, QString(), &ok).toStdString();
    if(ok){
        invokeOnConfig([this](ConfigState &c, const Plugin *plugin, const std::string &name, const Parameter &newParam){
            auto& map = c.parameters.at(plugin);
            if(map.find(name) == map.end()){
                // register new name in menu
                auto* fmtMenu = paramMenus.at(plugin);
                fmtMenu->addAction(name.c_str(), this, &MainWindow::loadParam);
            }
            // save parameter
            map[name] = newParam;
        }, curParam->getFmt(), name, *curParam);
    }
}

void MainWindow::savePreset()
{
    if(!presetWidget->curPreset){
        return;
    }
    const auto &curPreset = presetWidget->curPreset;
    bool ok;
    auto name = QInputDialog::getText(this, "Save IO preset", "Name of preset",
                                      QLineEdit::Normal, QString(), &ok).toStdString();
    if(ok){
        invokeOnConfig([this](ConfigState &c, const Plugin *plugin, const std::string &name, const Preset &newPreset){
            auto& map = c.presets.at(plugin);
            if(map.find(name) == map.end()){
                // register new name in menu
                auto* fmtMenu = presetMenus.at(plugin);
                fmtMenu->addAction(name.c_str(), this, &MainWindow::loadPreset);
            }
            // save preset
            map[name] = newPreset;
        }, curPreset->getFmt(), name, *curPreset);
    }
}

void MainWindow::saveScreenshot()
{
    QFileDialog diag{this};
    diag.setDirectory(path);
    diag.setAcceptMode(QFileDialog::AcceptSave);
    diag.setMimeTypeFilters({"image/png"});
    if(diag.exec() == QDialog::Accepted){
        path = diag.directory();
        auto target = diag.selectedFiles()[0];

        saveScreenshot(target);
    }
}

void MainWindow::saveScreenshots()
{
    QFileDialog diag{this};
    diag.setDirectory(path);
    diag.setFileMode(QFileDialog::Directory);
    diag.setAcceptMode(QFileDialog::AcceptSave);
    diag.setOption(QFileDialog::DontConfirmOverwrite);
    diag.setOption(QFileDialog::ShowDirsOnly, false);
    if(diag.exec() == QDialog::Accepted){
        path = diag.directory();
        size_t curStep = curVP->moldata[&curMol()].curStep;
        int width = std::log10(static_cast<double>(curMol().getNstep()))+1;
        for(size_t i=1; i <= curMol().getNstep(); ++i){
            curVP->setStep(i);
            saveScreenshot(fmt::format("{}/Step-{:0{}}.png", path.path().toStdString(), i, width).c_str());
        }
        curVP->setStep(curStep);
    }
}

void MainWindow::saveScreenshot(QString fn)
{
    if(!fn.endsWith(".png", Qt::CaseInsensitive)){
        fn += ".png";
    }
    // TODO: de-antialias for screenshot?
//    auto aa = _config.settings.antialias.val;
//    _config.settings.antialias.val = false;
    auto img = curVP->openGLWidget->grabFramebuffer();
    img.save(fn);
//    _config.settings.antialias.val = aa;
//    updateWidgets(0);
}

const ConfigState &MainWindow::config()
{
    return _config;
}
