#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent):
    QMainWindow(parent),
    ui(new Ui::MainWindow),
    curMol(NULL)
{
    ui->setupUi(this);
    connect(ui->actionAbout_Qt,SIGNAL(triggered()),qApp,SLOT(aboutQt()));
    QSignalBlocker tableBlocker(ui->cellVecTable);
    for(int j=0;j!=3;++j){
        for(int k=0;k!=3;++k){
            ui->cellVecTable->setItem(j,k,new QTableWidgetItem());
        }
    }
    newMol();
}

MainWindow::MainWindow(Vipster::Molecule m, QWidget *parent):
    QMainWindow(parent),
    ui(new Ui::MainWindow),
    curMol(NULL)
{
    ui->setupUi(this);
    connect(ui->actionAbout_Qt,SIGNAL(triggered()),qApp,SLOT(aboutQt()));
    QSignalBlocker tableBlocker(ui->cellVecTable);
    for(int j=0;j!=3;++j){
        for(int k=0;k!=3;++k){
             ui->cellVecTable->setItem(j,k,new QTableWidgetItem());
        }
    }
    newMol(m);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::setMol(int i)
{
    curMol = &molecules.at(i);
    uint steps = curMol->steps.size();
    uint curStep = curMol->stepIdx;
    //Step-control
    ui->stepLabel->setText(QString::number(steps));
    ui->stepEdit->setMaximum(steps);
    ui->stepEdit->setValue(curStep+1);
    ui->stepSlider->setMaximum(steps);
    ui->stepSlider->setValue(curStep+1);
    if(steps == 1){
        ui->stepEdit->setDisabled(true);
        ui->stepSlider->setDisabled(true);
    }else{
        ui->stepEdit->setEnabled(true);
        ui->stepSlider->setEnabled(true);
    }
    setStep(curStep);
}

void MainWindow::setStep(int i)
{
    if(i<0){
        i=curMol->stepIdx+1;
    }else{
        curMol->stepIdx=i;
    }
    //Handle control-buttons
    if(i == 1){
        ui->preStepButton->setDisabled(true);
        ui->firstStepButton->setDisabled(true);
    }else{
        ui->preStepButton->setEnabled(true);
        ui->firstStepButton->setEnabled(true);
    }
    //Fill atom list
    QSignalBlocker blockTable(ui->atomTable);
    int oldCount = ui->atomTable->rowCount();
    int nat = curMol->curStep().getNat();
    const std::vector<Vipster::Atom> *atoms = &curMol->curStep().getAtoms();
    ui->atomTable->setRowCount(nat);
    if( oldCount < nat){
        for(int j=oldCount;j!=nat;++j){
            for(int k=0;k!=4;++k){
                ui->atomTable->setItem(j,k,new QTableWidgetItem());
                ui->atomTable->item(j,k)->setFlags(
                            Qt::ItemIsSelectable|Qt::ItemIsEditable|
                            Qt::ItemIsUserCheckable|Qt::ItemIsEnabled);
            }
        }
    }
    for(int j=0;j!=nat;++j){
        const Vipster::Atom *at = &atoms->at(j);
        ui->atomTable->item(j,0)->setText(at->name.c_str());
        ui->atomTable->item(j,0)->setCheckState(Qt::CheckState(at->hidden*2));
        for(int k=0;k!=3;++k){
            ui->atomTable->item(j,k+1)->setText(QString::number(at->coord[k]));
            ui->atomTable->item(j,k+1)->setCheckState(Qt::CheckState(at->fix[k]*2));
        }
    }
    //Fill cell view
    ui->cellDimBox->setValue(curMol->curStep().getCellDim());
    std::array<Vipster::Vec,3> vec = curMol->curStep().getCellVec();
    for(int j=0;j!=3;++j){
        for(int k=0;k!=3;++k){
            ui->cellVecTable->item(j,k)->setText(QString::number(vec[j][k]));
        }
    }
    //Update display widget
    ui->openGLWidget->setStep(curMol->curStep());
}

void MainWindow::editAtoms()
{
    const QObject *sender = QObject::sender();
    if ( sender == ui->actionNew_Atom){
        curMol->curStep().newAtom();
    }
//    }else if ( sender == ui->actionDelete_Atom_s){
//        curMol->curStep().delAtom();
//    }
    setStep();
}

void MainWindow::newMol(Vipster::Molecule m)
{
    molecules.push_back(m);
    ui->molList->addItem(m.name.c_str());
    setMol(ui->molList->count()-1);
}

void MainWindow::about()
{
    QMessageBox::about(this,QString("About Vipster"),QString("aha?"));
}

void MainWindow::on_cellDimBox_valueChanged(double arg1)
{
    curMol->curStep().setCellDim(arg1,ui->cellScaleBox->isChecked());
    setStep();
}

void MainWindow::on_cellVecTable_cellChanged(int row, int column)
{
    std::array<std::array<float,3>,3> vec;
    vec = curMol->curStep().getCellVec();
    vec[row][column] = locale().toDouble(ui->cellVecTable->item(row,column)->text());
    curMol->curStep().setCellVec(vec);
    setStep();
}

void MainWindow::on_atomTable_cellChanged(int row, int column)
{
    Vipster::Atom at = curMol->curStep().getAtom(row);
    QTableWidgetItem *cell = ui->atomTable->item(row,column);
    switch(column){
    case 0:
        at.name = cell->text().toStdString();
        at.hidden = cell->checkState()/2;
        break;
    default:
        at.coord[column-1] = locale().toDouble(cell->text());
        at.fix[column-1] = cell->checkState()/2;
    }
    curMol->curStep().setAtom(row,at);
    setStep();
}
