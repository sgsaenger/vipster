#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QMessageBox>

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

MainWindow::MainWindow(const IO::BaseData &d, QWidget *parent):
    QMainWindow{parent},
    ui{new Ui::MainWindow}
{
    ui->setupUi(this);
    connect(ui->actionAbout_Qt,SIGNAL(triggered()),qApp,SLOT(aboutQt()));
    newMol(d.mol);
    newParam(d.param);
    setupUI();
}

MainWindow::MainWindow(IO::BaseData &&d, QWidget *parent):
    QMainWindow{parent},
    ui{new Ui::MainWindow}
{
    ui->setupUi(this);
    connect(ui->actionAbout_Qt,SIGNAL(triggered()),qApp,SLOT(aboutQt()));
    newMol(std::move(d.mol));
    newParam(std::move(d.param));
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
    tabifyDockWidget(ui->molDock, ui->paramDock);
    ui->molDock->raise();
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
    curParam = params.at(i).get();
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

void MainWindow::newMol(const Molecule &m)
{
    molecules.push_back(m);
    ui->molWidget->registerMol(m.getName());
}

void MainWindow::newMol(Molecule &&m)
{
    molecules.push_back(std::move(m));
    ui->molWidget->registerMol(molecules.back().getName());
}

void MainWindow::newParam(const std::unique_ptr<IO::BaseParam> &p)
{
    if(p){
        params.push_back(std::make_unique<IO::BaseParam>(*p));
        ui->paramWidget->registerParam(params.back()->name);
//        setParam(params.size());
    }
}

void MainWindow::newParam(std::unique_ptr<IO::BaseParam> &&p)
{
    if(p){
        params.push_back(std::move(p));
        ui->paramWidget->registerParam(params.back()->name);
//        setParam(params.size());
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
