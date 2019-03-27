#include "version.h"
#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "savefmtdialog.h"
#include "tools.h"

#include <QMessageBox>
#include <QFileDialog>
#include <QInputDialog>
#include <QApplication>
#include <QStyle>

using namespace Vipster;

MainWindow::MainWindow(QString path, QWidget *parent):
    QMainWindow{parent},
    ui{new Ui::MainWindow},
    path{path}
{
    ui->setupUi(this);
    connect(ui->actionAbout_Qt, &QAction::triggered, qApp, &QApplication::aboutQt);
    connect(&playTimer, &QTimer::timeout, ui->stepEdit, &QSpinBox::stepUp);
    setupUI();
    newMol();
}

MainWindow::MainWindow(QString path, std::vector<IO::Data> &&d, QWidget *parent):
    QMainWindow{parent},
    ui{new Ui::MainWindow},
    path{path}
{
    ui->setupUi(this);
    setupUI();
    connect(ui->actionAbout_Qt, &QAction::triggered, qApp, &QApplication::aboutQt);
    connect(&playTimer, &QTimer::timeout, ui->stepEdit, &QSpinBox::stepUp);
    for(auto&& mol: d){
        newData(std::move(mol));
    }
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::setupUI()
{
    ui->firstStepButton->setIcon(style()->standardIcon(QStyle::SP_MediaSkipBackward));
    ui->preStepButton->setIcon(style()->standardIcon(QStyle::SP_MediaSeekBackward));
    ui->playButton->setIcon(style()->standardIcon(QStyle::SP_MediaPlay));
    ui->nextStepButton->setIcon(style()->standardIcon(QStyle::SP_MediaSeekForward));
    ui->lastStepButton->setIcon(style()->standardIcon(QStyle::SP_MediaSkipForward));
#ifdef Q_OS_MACOS
    setDockOptions(dockOptions()^VerticalTabs);
#endif
    // setup left dock-area
    baseWidgets = {
        ui->paramDock,
        ui->configDock,
        ui->settingsDock,
        ui->mpseDock,
        ui->pseDock,
        ui->dataDock,
    };
    for(auto& w: baseWidgets){
        tabifyDockWidget(ui->molDock, w);
    }
    ui->molDock->raise();
    // setup right dock-area
    for(const auto& pair: makeToolWidgets(this)){
        auto *tmp = new QDockWidget(this);
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
        toolWidgets.push_back(tmp);
    }
    // fill in global PSE
    ui->pseWidget->setPSE(&Vipster::pse);
    // fill in menu-options
    for(const auto& pair:IOPlugins){
        const auto& param_map = Vipster::params[pair.first];
        if(!param_map.empty()){
            auto* param_menu = ui->menuLoad_Parameter_set->addMenu(
                        QString::fromStdString(pair.second->name));
            for(const auto& p: param_map){
                param_menu->addAction(QString::fromStdString(p.first),
                                      this, &MainWindow::loadParam);
            }
        }
        const auto& conf_map = Vipster::configs[pair.first];
        if(!conf_map.empty()){
            auto* conf_menu = ui->menuLoad_IO_Config->addMenu(
                        QString::fromStdString(pair.second->name));
            for(const auto& p: conf_map){
                conf_menu->addAction(QString::fromStdString(p.first),
                                     this, &MainWindow::loadConfig);
            }
        }
    }
}

void MainWindow::updateWidgets(guiChange_t change)
{
    // if necessary, make sure that data is up to date
    if(change & (GuiChange::atoms | GuiChange::fmt)){
        curStep->evaluateCache();
    }
    if(change & GuiChange::selection){
        curSel->evaluateCache();
    }
    // notify widgets
    ui->openGLWidget->updateWidget(change);
    ui->molWidget->updateWidget(change);
    for(auto& w: baseWidgets){
        static_cast<BaseWidget*>(w->widget())->updateWidget(change);
    }
    for(auto& w: toolWidgets){
        static_cast<BaseWidget*>(w->widget())->updateWidget(change);
    }
}

void MainWindow::setMol(int i)
{
    curMol = &*std::next(molecules.begin(), i);
    int nstep = static_cast<int>(curMol->getNstep());
    if(moldata[curMol].curStep < 0){
        moldata[curMol].curStep = nstep;
    }
    //Step-control
    ui->stepLabel->setText(QString::number(nstep));
    QSignalBlocker boxBlocker(ui->stepEdit);
    QSignalBlocker slideBlocker(ui->stepSlider);
    ui->stepEdit->setMaximum(nstep);
    ui->stepSlider->setMaximum(nstep);
    if(nstep == 1){
        ui->stepEdit->setDisabled(true);
        ui->stepSlider->setDisabled(true);
    }else{
        ui->stepEdit->setEnabled(true);
        ui->stepSlider->setEnabled(true);
    }
    setStep(moldata[curMol].curStep);
    updateWidgets(guiMolChanged);
}

void MainWindow::setStep(int i)
{
    moldata[curMol].curStep = i;
    curStep = &curMol->getStep(static_cast<size_t>(i-1));
    // if no previous selection exists, create one, afterwards assign it
    auto& tmpSel = stepdata[curStep].sel;
    if(!tmpSel && curStep){
        tmpSel = std::make_unique<Step::selection>(curStep->select(SelectionFilter{}));
    }
    curSel = tmpSel.get();
    //Handle control-elements
    if(playTimer.isActive() && (i == static_cast<int>(curMol->getNstep()))){
        playTimer.stop();
        ui->playButton->setIcon(style()->standardIcon(QStyle::SP_MediaPlay));
    }
    QSignalBlocker boxBlocker(ui->stepEdit);
    QSignalBlocker slideBlocker(ui->stepSlider);
    ui->stepEdit->setValue(i);
    ui->stepSlider->setValue(i);
    if(i == 1){
        ui->preStepButton->setDisabled(true);
        ui->firstStepButton->setDisabled(true);
    }else{
        ui->preStepButton->setEnabled(true);
        ui->firstStepButton->setEnabled(true);
    }
    if(i == static_cast<int>(curMol->getNstep())){
        ui->nextStepButton->setDisabled(true);
        ui->lastStepButton->setDisabled(true);
    }else{
        ui->nextStepButton->setEnabled(true);
        ui->lastStepButton->setEnabled(true);
    }
    //Update child widgets
    updateWidgets(guiStepChanged);
}

void MainWindow::stepBut(QAbstractButton* but)
{
    if(but == ui->firstStepButton){
        setStep(1);
    }else if(but == ui->lastStepButton){
        setStep(static_cast<int>(curMol->getNstep()));
    }else if(but == ui->playButton){
        if(playTimer.isActive()){
            playTimer.stop();
            ui->playButton->setIcon(style()->standardIcon(QStyle::SP_MediaPlay));
        }else{
            playTimer.start(static_cast<int>(settings.animstep.val));
            ui->playButton->setIcon(style()->standardIcon(QStyle::SP_MediaPause));
        }
    }
}

void MainWindow::editAtoms(QAction* sender)
{
    guiChange_t change{};
    if ( sender == ui->actionNew_Atom){
        curStep->newAtom();
        change = GuiChange::atoms;
    }else if ( sender == ui->actionDelete_Atom_s){
        curStep->delAtoms(*curSel);
        change = GuiChange::atoms | GuiChange::selection;
    }else if ( sender == ui->actionHide_Atom_s){
        for(auto& at: *curSel){
            at.properties->flags[AtomFlag::Hidden] = 1;
        }
        change = GuiChange::atoms;
    }else if ( sender == ui->actionShow_Atom_s){
        for(auto& at: *curSel){
            at.properties->flags[AtomFlag::Hidden] = 0;
        }
        change = GuiChange::atoms;
    }else if ( sender == ui->actionRename_Atom_s){
        auto tmp = QInputDialog::getText(this, "Rename atoms",
                                         "Enter new Atom-type for selected atoms:")
                   .toStdString();
        for(auto& at: *curSel){
            at.name = tmp;
        }
        change = GuiChange::atoms;
    }else if ( sender == ui->actionCopy_Atom_s){
        copyBuf = *curSel;
    }else if ( sender == ui->actionCut_Atom_s){
        copyBuf = *curSel;
        curStep->delAtoms(*curSel);
        change = GuiChange::atoms | GuiChange::selection;
    }else if ( sender == ui->actionPaste_Atom_s){
        curStep->newAtoms(copyBuf);
        change = GuiChange::atoms;
    }else if( sender == ui->actionSet_Bonds){
        curStep->setBonds();
        change = GuiChange::atoms;
    }
    if(change){
        updateWidgets(change);
    }
}

void MainWindow::newMol()
{
    molecules.emplace_back();
    ui->molWidget->registerMol(molecules.back().getName());
}

void MainWindow::newMol(QAction* sender)
{
    if(sender == ui->actionCopy_Trajector){
        molecules.emplace_back(*curMol);
        molecules.back().setName(molecules.back().getName() + " (copy)");
        ui->molWidget->registerMol(molecules.back().getName());
    }else if( sender == ui->actionCopy_single_Step){
        molecules.emplace_back(*curStep, curMol->getName() + " (copy of step " +
                               std::to_string(moldata[curMol].curStep) + ')');
        ui->molWidget->registerMol(molecules.back().getName());
    }else if( sender == ui->actionCopy_current_Selection){
        molecules.emplace_back(*curSel, curMol->getName() + " (copy of selection of step" +
                               std::to_string(moldata[curMol].curStep) + ')');
        ui->molWidget->registerMol(molecules.back().getName());
    }
}

void MainWindow::newData(IO::Data &&d)
{
    molecules.push_back(std::move(d.mol));
    ui->molWidget->registerMol(molecules.back().getName());
    if(d.param){
        ui->paramWidget->registerParam(std::move(d.param));
    }
    for(auto& dat: d.data){
        data.push_back(std::move(dat));
        ui->dataWidget->registerData(data.back()->name);
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
    std::map<std::string, int> ext2id;
    for(auto &iop: IOPlugins){
        ext2id[iop.second->extension] = formats.size();
        formats << QString::fromStdString(iop.second->name);
    }
    QString fmt_s;
    IOFmt fmt;
    bool got_fmt{false};
    if(fileDiag.exec() != 0){
        auto files = fileDiag.selectedFiles();
        path = fileDiag.directory();
        if(files.count() != 0){
            int current = 0;
            // try to determine filetype by extension
            auto first = files[0].toStdString();
            auto pos = first.find_last_of('.');
            if(pos != first.npos){
                auto ext = first.substr(pos+1);
                auto pos = ext2id.find(ext);
                if(pos != ext2id.end()){
                    current = pos->second;
                }
            }
            // request filetype
            fmt_s = QInputDialog::getItem(this, "Select format", "Format:",
                                          formats, current, false, &got_fmt);
            // load files
            if(got_fmt){
                fmt = static_cast<IOFmt>(formats.indexOf(fmt_s));
                for(auto &file: files){
                    IO::Data dat;
                    try {
                        dat = readFile(file.toStdString(), fmt);
                    } catch (const IO::Error& e) {
                        QMessageBox msg{this};
                        msg.setText(QString{"Could not open file \""}+file+"\":\n"+e.what());
                        msg.exec();
                        continue;
                    }
                    newData(std::move(dat));
                }
            }
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
                          IO::State{static_cast<size_t>(ui->stepSlider->value()-1),
                                    ui->molWidget->getAtomFmt(),
                                    ui->molWidget->getCellFmt()});
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
    return ui->paramWidget->params;
}

const decltype (ConfigWidget::configs)& MainWindow::getConfigs() const noexcept
{
    return ui->configWidget->configs;
}

void MainWindow::addExtraData(GUI::Data* dat)
{
    ui->openGLWidget->addExtraData(dat);
    updateWidgets(GuiChange::extra);
}

void MainWindow::delExtraData(GUI::Data* dat)
{
    ui->openGLWidget->delExtraData(dat);
    updateWidgets(GuiChange::extra);
}

const GUI::GlobalData& MainWindow::getGLGlobals()
{
    return ui->openGLWidget->globals;
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
    auto pos = Vipster::params[fmt].find(s->text().toStdString());
    if(pos != Vipster::params[fmt].end()){
        ui->paramWidget->registerParam(pos->second->copy());
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
    auto pos = Vipster::configs[fmt].find(s->text().toStdString());
    if(pos != Vipster::configs[fmt].end()){
        ui->configWidget->registerConfig(pos->second->copy());
    }else{
        throw Error("Invalid IO-config");
    }
}

void MainWindow::saveParam()
{
    if(!ui->paramWidget->curParam){
        return;
    }
    bool ok;
    auto name = QInputDialog::getText(this, "Save parameter set", "Name of preset",
                                      QLineEdit::Normal, QString(), &ok).toStdString();
    if(ok){
        Vipster::params[ui->paramWidget->curFmt][name] =
                ui->paramWidget->curParam->copy();
        Vipster::params[ui->paramWidget->curFmt][name]->name = name;
    }
}

void MainWindow::saveConfig()
{
    if(!ui->configWidget->curConfig){
        return;
    }
    bool ok;
    auto name = QInputDialog::getText(this, "Save IO-Config", "Name of preset",
                                      QLineEdit::Normal, QString(), &ok).toStdString();
    if(ok){
        Vipster::configs[ui->configWidget->curFmt][name] =
                ui->configWidget->curConfig->copy();
        Vipster::configs[ui->configWidget->curFmt][name]->name = name;
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
        auto img = ui->openGLWidget->grabFramebuffer();
        img.save(target);
        settings.antialias.val = aa;
        updateWidgets(0);
    }
}
