#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QMessageBox>

MainWindow::MainWindow(QWidget *parent):
    QMainWindow(parent),
    curMol(nullptr),
    curStep(nullptr),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    connect(ui->actionAbout_Qt,SIGNAL(triggered()),qApp,SLOT(aboutQt()));
    newMol();
}

MainWindow::MainWindow(Vipster::Molecule m, QWidget *parent):
    QMainWindow(parent),
    curMol(nullptr),
    curStep(nullptr),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    connect(ui->actionAbout_Qt,SIGNAL(triggered()),qApp,SLOT(aboutQt()));
    newMol(m);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::setMol(void)
{
    setMol(1);
}

void MainWindow::setMol(int i)
{
    curMol = &molecules.at(i-1);
    uint steps = curMol->getNstep();
    //Step-control
    ui->stepLabel->setText(QString::number(steps));
    ui->stepEdit->setMaximum(steps);
    ui->stepEdit->setValue(1);
    ui->stepSlider->setMaximum(steps);
    ui->stepSlider->setValue(1);
    if(steps == 1){
        ui->stepEdit->setDisabled(true);
        ui->stepSlider->setDisabled(true);
    }else{
        ui->stepEdit->setEnabled(true);
        ui->stepSlider->setEnabled(true);
    }
    setStep();
}

void MainWindow::setStep(void)
{
    setStep(1);
}

void MainWindow::setStep(int i)
{
    curStep = &curMol->getStep(i-1);
    //Handle control-buttons
    if(i == 1){
        ui->preStepButton->setDisabled(true);
        ui->firstStepButton->setDisabled(true);
    }else{
        ui->preStepButton->setEnabled(true);
        ui->firstStepButton->setEnabled(true);
    }
    //Update child widgets
    ui->openGLWidget->setStep(curStep);
    ui->molWidget->setStep(curStep);
}

void MainWindow::editAtoms()
{
    const QObject *sender = QObject::sender();
    if ( sender == ui->actionNew_Atom){
        curStep->newAtom();
    }
//    }else if ( sender == ui->actionDelete_Atom_s){
//        curMol->curStep().delAtom();
//    }
    setStep();
}

void MainWindow::newMol(Vipster::Molecule m)
{
    molecules.push_back(m);
    setMol(molecules.size());
}

void MainWindow::about()
{
    QMessageBox::about(this,QString("About Vipster"),QString("aha?"));
}
