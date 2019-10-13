#include "version.h"
#include "mainwindow.h"
#include "ui_mainwindow.h"
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

using namespace Vipster;

MainWindow::MainWindow(QString path, ConfigState& state,
                       std::vector<IO::Data> &&d, QWidget *parent):
    QMainWindow{parent},
    state{state},
    pte{std::get<0>(state)},
    settings{std::get<1>(state)},
    params{std::get<2>(state)},
    configs{std::get<3>(state)},
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
    hsplits.push_back(new QSplitter{vsplit});
    hsplits.front()->addWidget(viewports.front());
    setCentralWidget(vsplit);
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
        baseWidgets.push_back(pair.first);
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
        auto c = dynamic_cast<ConfigWidget*>(pair.first);
        if (c) configWidget = c;
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
    for(const auto& pair:IOPlugins){
        if(pair.second->arguments & IO::Plugin::Param){
            auto* param_menu = ui->menuLoad_Parameter_set->addMenu(
                        QString::fromStdString(pair.second->name));
            paramMenus[pair.first] = param_menu;
            const auto& param_map = params[pair.first];
            if(!param_map.empty()){
                for(const auto& p: param_map){
                    param_menu->addAction(QString::fromStdString(p.first),
                                          this, &MainWindow::loadParam);
                }
            }
        }
        if(pair.second->arguments & IO::Plugin::Config){
            auto* conf_menu = ui->menuLoad_IO_Config->addMenu(
                        QString::fromStdString(pair.second->name));
            configMenus[pair.first] = conf_menu;
            const auto& conf_map = configs[pair.first];
            if(!conf_map.empty()){
                for(const auto& p: conf_map){
                    conf_menu->addAction(QString::fromStdString(p.first),
                                         this, &MainWindow::loadConfig);
                }
            }
        }
    }
}

void MainWindow::updateWidgets(GUI::change_t change)
{
    // pull in mol/step selection from active viewport
    if((change & GUI::molChanged) == GUI::molChanged){
        curMol = viewports.front()->curMol;
    }
    if((change & GUI::stepChanged) == GUI::stepChanged){
        curStep = viewports.front()->curStep;
        curSel = viewports.front()->curSel;
    }
    // if necessary, make sure that data is up to date
    if(change & (GUI::Change::atoms | GUI::Change::fmt)){
        curStep->evaluateCache();
    }
    if(change & GUI::Change::selection){
        curSel->evaluateCache();
    }
    // notify widgets
    for(auto& w: viewports){
        w->updateWidget(change);
    }
    for(auto& w: baseWidgets){
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
            if(sender->active){
                throw Error{"Cannot close active viewport"};
            }else{
                auto pos = std::find(viewports.begin(), viewports.end(), sender);
                if(pos == viewports.end()){
                    throw Error{"Invalid viewport"};
                }else{
                    viewports.erase(pos);
                    delete sender;
                }
            }
            break;
        case VP_VSPLIT:
            viewports.push_back(new ViewPort{this});
            viewports.back()->updateWidget(GUI::stepChanged|GUI::molChanged);
            hsplits.push_back(new QSplitter{vsplit});
            hsplits.back()->addWidget(viewports.back());
            vsplit->addWidget(hsplits.back());
            break;
        case VP_HSPLIT:
            viewports.push_back(new ViewPort{this});
            viewports.back()->updateWidget(GUI::stepChanged|GUI::molChanged);
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
            break;
        default:
            throw Error{"Invalid viewport change"};
        }
    }catch(const Error& e){
        QMessageBox::information(this, "ViewPort Error", e.what());
    }
}

void MainWindow::setBondMode(bool b)
{
    // TODO: keep tied viewports in sync?
    viewports.front()->setBondMode(b);
}

void MainWindow::setMultEnabled(bool b)
{
    // TODO: keep tied viewports in sync?
    viewports.front()->setMultEnabled(b);
}

void MainWindow::editAtoms(QAction* sender)
{
    GUI::change_t change{};
    if ( sender == ui->actionNew_Atom){
        curStep->newAtom();
        change = GUI::Change::atoms;
    }else if ( sender == ui->actionDelete_Atom_s){
        curStep->delAtoms(*curSel);
        change = GUI::Change::atoms | GUI::Change::selection;
    }else if ( sender == ui->actionHide_Atom_s){
        for(auto& at: *curSel){
            at.properties->flags[AtomFlag::Hidden] = 1;
        }
        change = GUI::Change::atoms;
    }else if ( sender == ui->actionShow_Atom_s){
        for(auto& at: *curSel){
            at.properties->flags[AtomFlag::Hidden] = 0;
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
        copyBuf = *curSel;
    }else if ( sender == ui->actionCut_Atom_s){
        copyBuf = *curSel;
        curStep->delAtoms(*curSel);
        change = GUI::Change::atoms | GUI::Change::selection;
    }else if ( sender == ui->actionPaste_Atom_s){
        curStep->newAtoms(copyBuf);
        change = GUI::Change::atoms;
    }
    if(change){
        updateWidgets(change);
    }
}

void MainWindow::registerMol(const std::string& name)
{
    // delegate to active viewport which has relevant molList
    for(auto& w: viewports){
        w->registerMol(name);
    }
}

void MainWindow::newMol()
{
    molecules.emplace_back();
    molecules.back().pte->root = &pte;
    registerMol(molecules.back().getName());
}

void MainWindow::newMol(QAction* sender)
{
    if(sender == ui->actionCopy_Trajector){
        molecules.emplace_back(*curMol);
        molecules.back().setName(molecules.back().getName() + " (copy)");
        molecules.back().pte->root = &pte;
        registerMol(molecules.back().getName());
    }else if( sender == ui->actionCopy_single_Step){
        molecules.emplace_back(*curStep, curMol->getName() + " (copy of step " +
                               std::to_string(moldata[curMol].curStep) + ')');
        molecules.back().pte->root = &pte;
        registerMol(molecules.back().getName());
    }else if( sender == ui->actionCopy_current_Selection){
        molecules.emplace_back(*curSel, curMol->getName() + " (copy of selection of step" +
                               std::to_string(moldata[curMol].curStep) + ')');
        molecules.back().pte->root = &pte;
        registerMol(molecules.back().getName());
    }
}

void MainWindow::newData(IO::Data &&d)
{
    molecules.push_back(std::move(d.mol));
    molecules.back().pte->root = &pte;
    registerMol(molecules.back().getName());
    if(d.param){
        paramWidget->registerParam(std::move(d.param));
    }
    for(auto& dat: d.data){
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
    for(auto &iop: IOPlugins){
        formats << QString::fromStdString(iop.second->name);
    }
    if(fileDiag.exec() != 0){
        auto files = fileDiag.selectedFiles();
        path = fileDiag.directory();
        // if cancelled by user, just return
        if(files.empty()) return;
        // else try to read file
        const auto& file = files[0].toStdString();
        try {
            try {
                // try to open file without explicit format
                newData(readFile(file));
            } catch (const IO::Error &e) {
                // if error is not fatal, we just need the format, so request it
                if(!e.fatal){
                    bool got_fmt{false};
                    auto fmt_s = QInputDialog::getItem(this, "Select format", "Format:",
                                                       formats, 0, false, &got_fmt);
                    // if the user selected the format, read the file
                    if(got_fmt){
                        auto fmt = static_cast<IOFmt>(formats.indexOf(fmt_s));
                        newData(readFile(file, fmt));
                    }
                }else{
                    throw;
                }
            }
        }catch (const IO::Error &e){
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
        SaveFmtDialog sfd{this};
        if(sfd.exec() == QDialog::Accepted){
            try{
                writeFile(target, sfd.fmt, *curMol,
                          sfd.getParam(), sfd.getConfig(),
                          moldata[curMol].curStep);
            }catch(const IO::Error& e){
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

const decltype (ConfigWidget::configs)& MainWindow::getConfigs() const noexcept
{
    return configWidget->configs;
}

void MainWindow::addExtraData(GUI::Data* dat)
{
//    ui->openGLWidget->addExtraData(dat);
    updateWidgets(GUI::Change::extra);
}

void MainWindow::delExtraData(GUI::Data* dat)
{
//    ui->openGLWidget->delExtraData(dat);
    updateWidgets(GUI::Change::extra);
}

const GUI::GlobalData& MainWindow::getGLGlobals()
{
//    return ui->openGLWidget->globals;
    return GUI::GlobalData{};
}

void MainWindow::loadParam()
{
    auto* s = static_cast<QAction*>(sender());
    auto* p = static_cast<QMenu*>(s->parent());
    IOFmt fmt = [](QString name){
        for(const auto& pair:IOPlugins){
            if(name == pair.second->name.c_str()){
                return pair.first;
            }
        }
        throw Error("Invalid parameter set");
    }(p->title());
    auto pos = params[fmt].find(s->text().toStdString());
    if(pos != params[fmt].end()){
        paramWidget->registerParam(pos->second->copy());
    }else{
        throw Error("Invalid parameter set");
    }
}

void MainWindow::loadConfig()
{
    auto* s = static_cast<QAction*>(sender());
    auto* p = static_cast<QMenu*>(s->parent());
    IOFmt fmt = [](QString name){
        for(const auto& pair:IOPlugins){
            if(name == pair.second->name.c_str()){
                return pair.first;
            }
        }
        throw Error("Invalid IO-config");
    }(p->title());
    auto pos = configs[fmt].find(s->text().toStdString());
    if(pos != configs[fmt].end()){
        configWidget->registerConfig(pos->second->copy());
    }else{
        throw Error("Invalid IO-config");
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
        map[name] = curParam->copy();
        map[name]->name = name;
    }
}

void MainWindow::saveConfig()
{
    if(!configWidget->curConfig){
        return;
    }
    const auto &curConfig = configWidget->curConfig;
    bool ok;
    auto name = QInputDialog::getText(this, "Save IO-Config", "Name of preset",
                                      QLineEdit::Normal, QString(), &ok).toStdString();
    if(ok){
        auto& map = configs[curConfig->getFmt()];
        if(map.find(name) == map.end()){
            // register new name in menu
            auto* fmtMenu = configMenus.at(curConfig->getFmt());
            fmtMenu->addAction(name.c_str(), this, &MainWindow::loadConfig);
        }
        // save config
        map[name] = curConfig->copy();
        map[name]->name = name;
    }
}

void MainWindow::about()
{
    QMessageBox::about(this,QString("About Vipster"),
    QString("<h2>Vipster v" VIPSTER_VERSION "b</h2>"
            "<p>"
            "©Sebastian Gsänger, 2019"
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
        auto target = diag.selectedFiles()[0];
        if(!target.endsWith(".png", Qt::CaseInsensitive)){
            target += ".png";
        }
        path = diag.directory();

        auto aa = settings.antialias.val;
        settings.antialias.val = false;
        auto img = viewports.front()->openGLWidget->grabFramebuffer();
        img.save(target);
        settings.antialias.val = aa;
        updateWidgets(0);
    }
}
