#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "savefmtdialog.h"
#include <QMessageBox>
#include <QFileDialog>
#include <QInputDialog>

using namespace Vipster;

MainWindow::MainWindow(QWidget *parent):
    QMainWindow{parent},
    ui{new Ui::MainWindow}
{
    ui->setupUi(this);
    connect(ui->actionAbout_Qt,SIGNAL(triggered()),qApp,SLOT(aboutQt()));
    newMol();
    setupUI();
}

MainWindow::MainWindow(std::vector<IO::Data> &&d, QWidget *parent):
    QMainWindow{parent},
    ui{new Ui::MainWindow}
{
    ui->setupUi(this);
    connect(ui->actionAbout_Qt,SIGNAL(triggered()),qApp,SLOT(aboutQt()));
    for(auto&& mol: d){
        newData(std::move(mol));
    }
    setupUI();
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
    // setup left dock-area
    tabifyDockWidget(ui->molDock, ui->paramDock);
    ui->molDock->raise();
#ifdef Q_OS_MACOS
    setDockOptions(dockOptions()^VerticalTabs);
#endif
    // fill in menu-options
    ui->menuPWScf->removeAction(ui->actionParamDummy); // qtcreator only creates menu when actions are provided
    auto pwi_range = Vipster::params.equal_range(IOFmt::PWI);
    for(auto it=pwi_range.first; it!=pwi_range.second; ++it){
        ui->menuPWScf->addAction(QString::fromStdString(it->second->name),
                                 this, &MainWindow::loadParam);
    }
}

void MainWindow::updateWidgets(Change change)
{
    if(change & (Change::atoms | Change::cell))
        ui->openGLWidget->setStep(curStep);
    ui->molWidget->updateWidget(change);
    ui->paramWidget->updateWidget(change);
}

void MainWindow::setFmt(int i, bool apply, bool scale)
{
    fmt = static_cast<AtomFmt>(i);
    if(apply){
        curStep->setFmt(fmt, scale);
        updateWidgets(stepChanged);
    }else{
        updateWidgets(Change::fmt);
    }
}

AtomFmt MainWindow::getFmt()
{
    return fmt;
}

void MainWindow::setMol(int i)
{
    curMol = &molecules.at(i);
    uint steps = curMol->getNstep();
    //Step-control
    ui->stepLabel->setText(QString::number(steps));
    QSignalBlocker boxBlocker(ui->stepEdit);
    QSignalBlocker slideBlocker(ui->stepSlider);
    ui->stepEdit->setMaximum(steps);
    ui->stepEdit->setValue(steps);
    ui->stepSlider->setMaximum(steps);
    ui->stepSlider->setValue(steps);
    if(steps == 1){
        ui->stepEdit->setDisabled(true);
        ui->stepSlider->setDisabled(true);
    }else{
        ui->stepEdit->setEnabled(true);
        ui->stepSlider->setEnabled(true);
    }
    setStep(steps);
    updateWidgets(molChanged);
}

void MainWindow::setStep(int i)
{
    curStep = &curMol->getStep(i-1);
    fmt = curStep->getFmt();
    //Handle control-buttons
    if(i == 1){
        ui->preStepButton->setDisabled(true);
        ui->firstStepButton->setDisabled(true);
    }else{
        ui->preStepButton->setEnabled(true);
        ui->firstStepButton->setEnabled(true);
    }
    if(i == curMol->getNstep()){
        ui->nextStepButton->setDisabled(true);
        ui->lastStepButton->setDisabled(true);
    }else{
        ui->nextStepButton->setEnabled(true);
        ui->lastStepButton->setEnabled(true);
    }
    //Update child widgets
    updateWidgets(stepChanged);
}

void MainWindow::setParam(int i)
{
    curParam = params.at(i).second.get();
    updateWidgets(Change::param);
}

void MainWindow::stepBut(QAbstractButton* but)
{
    if(but == ui->firstStepButton){
        setStep(1);
    }else if(but == ui->lastStepButton){
        setStep(curMol->getNstep());
    }else if(but == ui->playButton){
        //TODO
    }
}

void MainWindow::editAtoms()
{
    const QObject *sender = QObject::sender();
    if ( sender == ui->actionNew_Atom){
        curStep->newAtom();
    }
    //TODO:
//    }else if ( sender == ui->actionDelete_Atom_s){
//        curMol->curStep().delAtom();
//    }
    updateWidgets(Change::atoms);
}

void MainWindow::newMol()
{
    molecules.emplace_back();
    ui->molWidget->registerMol(molecules.back().getName());
}

void MainWindow::newData(IO::Data &&d)
{
//    newMol(std::move(d.mol));
    molecules.push_back(std::move(d.mol));
    ui->molWidget->registerMol(molecules.back().getName());
//    newParam(d.fmt, std::move(d.param));
    if(d.param){
        params.push_back({d.fmt, std::move(d.param)});
        ui->paramWidget->registerParam(d.fmt, params.back().second->name);
    }
}

void MainWindow::loadMol()
{
    // File dialog
    QStringList files;
    QFileDialog fileDiag{this};
    fileDiag.setDirectory(path);
    fileDiag.setFileMode(QFileDialog::ExistingFiles);
    // Format dialog
    QStringList formats{};
    for(auto &iop: IOPlugins){
        formats << QString::fromStdString(iop.second->name);
    }
    QString fmt_s;
    IOFmt fmt;
    bool got_fmt{false};
    if(fileDiag.exec()){
        files = fileDiag.selectedFiles();
        path = fileDiag.directory();
        if(files.count()){
            fmt_s = QInputDialog::getItem(this, "Select format", "Format:",
                                          formats, 0, false, &got_fmt);
            if(got_fmt){
                fmt = static_cast<IOFmt>(formats.indexOf(fmt_s));
                for(auto &file: files){
                    newData(readFile(file.toStdString(), fmt));
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
    if(fileDiag.exec()){
        auto target = fileDiag.selectedFiles()[0].toStdString();
        path = fileDiag.directory();
        SaveFmtDialog sfd{this};
        if(sfd.exec()){
            writeFile(target, sfd.fmt, *curMol, sfd.param, sfd.config);
        }
    }
}

void MainWindow::loadParam()
{
    auto* s = static_cast<QAction*>(sender());
    auto* p = static_cast<QMenu*>(s->parent());
    IOFmt fmt = [&](){
        if(p->title() == "PWScf"){
            return IOFmt::PWI;
        }
        throw Error("Invalid parameter set");
    }();
    auto range = Vipster::params.equal_range(fmt);
    for(auto it=range.first; it!=range.second; ++it){
        if (it->second->name.c_str() == s->text()){
            params.push_back({fmt, it->second->copy()});
            ui->paramWidget->registerParam(fmt, params.back().second->name);
            return;
        }
    }
}

void MainWindow::about()
{
    QMessageBox::about(this,QString("About Vipster"),
    QString("<h2>Vipster</h2>"
            "<p>"
            "©Sebastian Gsänger, 2018"
            "<br>"
            "<a href='https://hein09.github.io/vipster'>Homepage</a>"
            "<br>"
            "<a href='https://github.com/hein09/vipster'>Source</a>"
            "</p>"
            "<p>"
            "This program is provided under the GPLv3."
            "<br>"
            "It uses"
            "<a href='https://github.com/nlohmann/json'>JSON for Modern C++</a>, "
            "<a href='https://github.com/catchorg/catch2'>Catch2</a> and "
            "<a href='https://github.com/pybind/pybind11'>pybind11</a>."
            "</p>"));
}
