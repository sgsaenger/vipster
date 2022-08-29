#include "version.h"
#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "ui_viewport.h"
#include "savefmtdialog.h"
#include "mainwidgets.h"
#include "toolwidgets.h"
#include "glwidget.h"

#include <QDockWidget>
#include <QMessageBox>
#include <QFileDialog>
#include <QInputDialog>
#include <QApplication>
#include <QStyle>

#include <filesystem>
#include <fmt/format.h>

using namespace Vipster;
namespace fs = std::filesystem;

MainWindow::MainWindow(QString path, ConfigState& state,
                       std::vector<IOTuple> &&d, QWidget *parent):
    QMainWindow{parent},
    state{state},
    pte{std::get<0>(state)},
    settings{std::get<1>(state)},
    plugins{std::get<2>(state)},
    params{std::get<3>(state)},
    presets{std::get<4>(state)},
    ui{new Ui::MainWindow},
    path{path}
{
    ui->setupUi(this);
    setupUI();
    connect(ui->actionAbout_Qt, &QAction::triggered, qApp, &QApplication::aboutQt);
    // create first level of splitters and main viewport
    viewports.push_back(new ViewPort{this, true});
    vsplit = new QSplitter{this};
    vsplit->setOrientation(Qt::Vertical);
    vsplit->setChildrenCollapsible(false);
    hsplits.push_back(new QSplitter{vsplit});
    hsplits.back()->setChildrenCollapsible(false);
    hsplits.back()->addWidget(viewports.front());
    setCentralWidget(vsplit);
    curVP = viewports.front();
    curVP->makeActive(true);
    curVP->ui->closeButton->setDisabled(true);
    // load molecules
    if(d.empty()){
        newMol();
    }else{
        for(auto&& mol: d){
            newData(std::move(mol));
        }
    }
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::setupUI()
{
#ifdef Q_OS_MACOS
    setDockOptions(dockOptions()^VerticalTabs);
#endif
    // setup left dock-area
    QDockWidget* firstDock{nullptr};
    for(const auto& pair: makeMainWidgets(this)){
        auto *tmp = new QDockWidget(this);
        mainWidgets.push_back(pair.first);
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
    // setup right dock-area
    for(const auto& pair: makeToolWidgets(this)){
        auto *tmp = new QDockWidget(this);
        toolWidgets.push_back(pair.first);
        tmp->setWidget(pair.first);
        tmp->setAllowedAreas(Qt::BottomDockWidgetArea|
                             Qt::RightDockWidgetArea|
                             Qt::TopDockWidgetArea);
        tmp->setWindowTitle(pair.second);
        addDockWidget(Qt::RightDockWidgetArea, tmp);
        auto action = ui->toolBar->addAction(pair.second);
        action->setCheckable(true);
        connect(action, &QAction::toggled, tmp, &QWidget::setVisible);
        tmp->hide();
    }
    // fill in menu-options
    for(const auto& plug: plugins){
        if(plug->makeParam){
            auto* param_menu = ui->menuLoad_Parameter_set->addMenu(
                        QString::fromStdString(plug->name));
            paramMenus[plug] = param_menu;
            const auto& param_map = params[plug];
            if(!param_map.empty()){
                for(const auto& p: param_map){
                    param_menu->addAction(QString::fromStdString(p.first),
                                          this, &MainWindow::loadParam);
                }
            }
        }
        if(plug->makePreset){
            auto* conf_menu = ui->menuLoad_IO_Preset->addMenu(
                        QString::fromStdString(plug->name));
            presetMenus[plug] = conf_menu;
            const auto& conf_map = presets[plug];
            if(!conf_map.empty()){
                for(const auto& p: conf_map){
                    conf_menu->addAction(QString::fromStdString(p.first),
                                         this, &MainWindow::loadPreset);
                }
            }
        }
    }
}

void MainWindow::updateWidgets(GUI::change_t change)
{
    // pull in mol/step selection from active viewport
    if((change & GUI::molChanged) == GUI::molChanged){
        curMol = curVP->curMol;
    }
    if((change & GUI::stepChanged) == GUI::stepChanged){
        curStep = curVP->curStep;
        curSel = curVP->curSel;
    }
    // if necessary, make sure that bonds/overlaps are up to date
    if((change & GUI::Change::atoms) &&
       (stepdata[curStep].automatic_bonds ||
        settings.overlap.val)){
        curStep->generateBonds(!stepdata[curStep].automatic_bonds);
    }
    // notify widgets
    for(auto& w: viewports){
        w->updateWidget(change);
    }
    for(auto& w: mainWidgets){
        w->updateWidget(change);
    }
    for(auto& w: toolWidgets){
        w->updateWidget(change);
    }
}

void MainWindow::changeViewports(ViewPort *sender, VPChange change)
{
    try{
        switch(change){
        case VP_CLOSE:
            {
                bool wasActive = sender->active;
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
            break;
        case VP_VSPLIT:
            if(viewports.size() == 1) viewports.front()->ui->closeButton->setEnabled(true);
            viewports.push_back(new ViewPort{*sender});
            viewports.back()->updateWidget(GUI::molChanged);
            hsplits.push_back(new QSplitter{vsplit});
            hsplits.back()->setChildrenCollapsible(false);
            hsplits.back()->addWidget(viewports.back());
            vsplit->addWidget(hsplits.back());
            break;
        case VP_HSPLIT:
            if(viewports.size() == 1) viewports.front()->ui->closeButton->setEnabled(true);
            viewports.push_back(new ViewPort{*sender});
            viewports.back()->updateWidget(GUI::molChanged);
            for(auto& h: hsplits){
                for(auto& c: h->children()){
                    if (c == sender){
                        h->addWidget(viewports.back());
                        return;
                    }
                }
            }
            throw Error{"Could not determine horizontal splitter for viewport"};
        case VP_ACTIVE:
            // deactivate other viewports
            for(auto& vp: viewports){
                vp->makeActive(vp == sender);
            }
            // make sender current viewport
            curVP = sender;
            updateWidgets(GUI::molChanged);
            break;
        default:
            throw Error{"Invalid viewport change"};
        }
    }catch(const Error& e){
        QMessageBox::information(this, "ViewPort Error", e.what());
    }
}

void MainWindow::editAtoms(QAction* sender)
{
    GUI::change_t change{};
    if ( sender == ui->actionNew_Atom){
        curStep->newAtom("C");
        change = GUI::Change::atoms;
    }else if ( sender == ui->actionDelete_Atom_s){
        curStep->delAtoms(*curSel);
        *curSel = curStep->select({});
        change = GUI::Change::atoms | GUI::Change::selection;
    }else if ( sender == ui->actionHide_Atom_s){
        for(auto& at: *curSel){
            at.properties->flags[AtomProperties::Hidden] = 1;
        }
        change = GUI::Change::atoms;
    }else if ( sender == ui->actionShow_Atom_s){
        for(auto& at: *curSel){
            at.properties->flags[AtomProperties::Hidden] = 0;
        }
        change = GUI::Change::atoms;
    }else if ( sender == ui->actionRename_Atom_s){
        auto tmp = QInputDialog::getText(this, "Rename atoms",
                                         "Enter new Atom-type for selected atoms:")
                   .toStdString();
        for(auto& at: *curSel){
            at.name = tmp;
        }
        change = GUI::Change::atoms;
    }else if ( sender == ui->actionCopy_Atom_s){
        copyBuf = std::make_unique<Step::selection>(*curSel);
    }else if ( sender == ui->actionCut_Atom_s){
        copyBuf = std::make_unique<Step::selection>(*curSel);
        curStep->delAtoms(*curSel);
        *curSel = curStep->select({});
        change = GUI::Change::atoms | GUI::Change::selection;
    }else if ( sender == ui->actionPaste_Atom_s){
        if(!copyBuf) return;
        auto oldNat = curStep->getNat();
        curStep->newAtoms(*copyBuf);
        *curVP->stepdata[curStep].sel = curStep->select(
            "index "+std::to_string(oldNat)+'-'+std::to_string(curStep->getNat()-1));
        change = GUI::Change::atoms | GUI::Change::selection;
    }
    if(change){
        updateWidgets(change);
    }
}

void MainWindow::registerMol(const std::string& name)
{
    for(auto& w: viewports){
        w->registerMol(name);
    }
}

void MainWindow::newMol()
{
    newMol(Molecule{});
}

void MainWindow::newMol(Molecule&& mol)
{
    molecules.emplace_back(std::move(mol));
    molecules.back().getPTE().root = &pte;
    registerMol(molecules.back().name);
}

void MainWindow::newMol(QAction* sender)
{
    if(sender == ui->actionCopy_Trajector){
        molecules.emplace_back(*curMol);
        molecules.back().name += " (copy)";
        molecules.back().getPTE().root = &pte;
        registerMol(molecules.back().name);
    }else if( sender == ui->actionCopy_single_Step){
        molecules.emplace_back(*curStep, curMol->name + " (copy of step " +
                               std::to_string(curVP->moldata[curMol].curStep) + ')');
        molecules.back().getPTE().root = &pte;
        registerMol(molecules.back().name);
    }else if( sender == ui->actionCopy_current_Selection){
        molecules.emplace_back(*curSel, curMol->name + " (copy of selection of step " +
                               std::to_string(curVP->moldata[curMol].curStep) + ')');
        molecules.back().getPTE().root = &pte;
        registerMol(molecules.back().name);
    }
}

void MainWindow::newData(IOTuple &&d)
{
    newMol(std::move(std::get<0>(d)));
    const auto& name = molecules.back().name;
    if(auto &param = std::get<1>(d)){
        paramWidget->registerParam(name, std::move(*param));
    }
    for(auto& dat: std::get<2>(d)){
        data.push_back(std::move(dat));
        updateWidgets(GUI::Change::data);
    }
}

void MainWindow::loadMol()
{
    // File dialog
    QFileDialog fileDiag{this};
    fileDiag.setDirectory(path);
    // TODO: limited to one file for now
    fileDiag.setFileMode(QFileDialog::ExistingFile);
    // Format dialog
    QStringList formats{};
    for(const auto plug: plugins){
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
        auto plugin = guessFmt(file, std::get<2>(state));
        if (!plugin){
            bool got_fmt{false};
            auto fmt_s = QInputDialog::getItem(this, "Select format", "Format:",
                                               formats, 0, false, &got_fmt);
            // if the user selected the format, read the file
            if(got_fmt){
                auto fmt = std::find_if(plugins.begin(), plugins.end(),
                    [&](const auto& plug){return plug->name.c_str() == fmt_s;});
                if(fmt == plugins.end())
                    throw Error{"Invalid format in loadMol occured"};
                plugin = *fmt;
            }
        }
        try {
            if(plugin){
                // try to open file
                newData(readFile(file, plugin));
            }
        }catch (const IOError &e){
            QMessageBox msg{this};
            msg.setText(QString{"Could not open file \""}+file.c_str()+"\":\n"+e.what());
            msg.exec();
        }
    }
}

void MainWindow::saveMol()
{
    QFileDialog fileDiag{this};
    fileDiag.setDirectory(path);
    fileDiag.setAcceptMode(QFileDialog::AcceptSave);
    if(fileDiag.exec() == QDialog::Accepted){
        auto target = fileDiag.selectedFiles()[0].toStdString();
        path = fileDiag.directory();
        SaveFmtDialog sfd{plugins, this};
        if(sfd.exec() == QDialog::Accepted){
            try{
                writeFile(target, sfd.plugin, *curMol,
                          curVP->moldata[curMol].curStep-1,
                          sfd.getParam(), sfd.getPreset());
            }catch(const IOError& e){
                QMessageBox msg{this};
                msg.setText(QString{"Could not write file \""}+target.c_str()+"\":\n"+e.what());
                msg.exec();
            }
        }
    }
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
    auto fmt = std::find_if(plugins.begin(), plugins.end(),
                 [&](const auto& plug){return plug->name.c_str() == p->title();});
    if(fmt == plugins.end()){
        throw Error{"Invalid parameter set"};
    }
    auto name = s->text().toStdString();
    auto pos = params[*fmt].find(name);
    if(pos != params[*fmt].end()){
        paramWidget->registerParam(name, pos->second);
    }else{
        throw Error("Invalid parameter set");
    }
}

void MainWindow::loadPreset()
{
    auto* s = static_cast<QAction*>(sender());
    auto* p = static_cast<QMenu*>(s->parent());
    auto fmt = std::find_if(plugins.begin(), plugins.end(),
                 [&](const auto& plug){return plug->name.c_str() == p->title();});
    if(fmt == plugins.end()){
        throw Error{"Invalid IO-preset"};
    }
    auto name = s->text().toStdString();
    auto pos = presets[*fmt].find(name);
    if(pos != presets[*fmt].end()){
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
        auto& map = params[curParam->getFmt()];
        if(map.find(name) == map.end()){
            // register new name in menu
            auto* fmtMenu = paramMenus.at(curParam->getFmt());
            fmtMenu->addAction(name.c_str(), this, &MainWindow::loadParam);
        }
        // save parameter
        map[name] = *curParam;
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
        auto& map = presets[curPreset->getFmt()];
        if(map.find(name) == map.end()){
            // register new name in menu
            auto* fmtMenu = presetMenus.at(curPreset->getFmt());
            fmtMenu->addAction(name.c_str(), this, &MainWindow::loadPreset);
        }
        // save preset
        map[name] = *curPreset;
    }
}

void MainWindow::about()
{
    QMessageBox::about(this,QString("About Vipster"),
    QString("<h2>Vipster v" VIPSTER_VERSION " (" VIPSTER_PLATFORM " " VIPSTER_ARCH ")</h2>"
            "<p>"
            "©Sebastian Gsänger, 2022"
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
        size_t curStep = curVP->moldata[curMol].curStep;
        int width = std::log10(static_cast<double>(curMol->getNstep()))+1;
        for(size_t i=1; i<=curMol->getNstep(); ++i){
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
    auto aa = settings.antialias.val;
    settings.antialias.val = false;
    auto img = curVP->openGLWidget->grabFramebuffer();
    img.save(fn);
    settings.antialias.val = aa;
    updateWidgets(0);
}
