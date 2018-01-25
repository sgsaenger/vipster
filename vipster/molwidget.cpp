#include "molwidget.h"
#include "ui_molwidget.h"
#include "atom.h"
#include <QTableWidgetItem>

MolWidget::MolWidget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::MolWidget)
{
    ui->setupUi(this);
    QSignalBlocker tableBlocker(ui->cellVecTable);
    for(int j=0;j!=3;++j){
        for(int k=0;k!=3;++k){
             ui->cellVecTable->setItem(j,k,new QTableWidgetItem());
        }
    }
}

MolWidget::~MolWidget()
{
    delete ui;
}

void MolWidget::setStep(Vipster::Step *step)
{
    curStep = step;
    //Fill atom list
    QSignalBlocker blockTable(ui->atomTable);
    QSignalBlocker blockCell(ui->cellVecTable);
    QSignalBlocker blockDim(ui->cellDimBox);
    QSignalBlocker blockEnabled(ui->cellEnabled);
    int oldCount = ui->atomTable->rowCount();
    int nat = curStep->getNat();
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
        const Vipster::Atom at = (*curStep)[j];
        //TODO
        ui->atomTable->item(j,0)->setText(at.name.c_str());
        ui->atomTable->item(j,0)->setCheckState(Qt::CheckState(at.hidden*2));
        for(int k=0;k!=3;++k){
            ui->atomTable->item(j,k+1)->setText(QString::number(at.coord[k]));
            ui->atomTable->item(j,k+1)->setCheckState(Qt::CheckState(at.fix[k]*2));
        }
    }
    //Fill cell view
    ui->cellEnabled->setChecked(curStep->hasCell());
    ui->cellDimBox->setValue( curStep->getCellDim(
            (Vipster::CdmFmt)ui->cellFmt->currentIndex()));
    Vipster::Mat vec = curStep->getCellVec();
    for(int j=0;j!=3;++j){
        for(int k=0;k!=3;++k){
            ui->cellVecTable->item(j,k)->setText(QString::number(vec[j][k]));
        }
    }
}

void MolWidget::on_cellEnabled_toggled(bool checked)
{
    curStep->enableCell(checked);
    emit stepChanged();
}

void MolWidget::on_cellFmt_currentIndexChanged(int idx)
{
    QSignalBlocker blockCDB(ui->cellDimBox);
    ui->cellDimBox->setValue(curStep->getCellDim((Vipster::CdmFmt)idx));
}

void MolWidget::on_cellDimBox_valueChanged(double cdm)
{
    curStep->setCellDim(cdm, (Vipster::CdmFmt)ui->cellFmt->currentIndex(), ui->cellScaleBox->isChecked());
    emit stepChanged();
}

void MolWidget::on_cellVecTable_cellChanged(int row, int column)
{
    Vipster::Mat vec;
    vec = curStep->getCellVec();
    vec[row][column] = locale().toDouble(ui->cellVecTable->item(row,column)->text());
    curStep->setCellVec(vec, ui->cellScaleBox->isChecked());
    emit stepChanged();
}

void MolWidget::on_atomTable_cellChanged(int row, int column)
{
    Vipster::Atom at = (*curStep)[row];
    const QTableWidgetItem *cell = ui->atomTable->item(row,column);
    switch(column){
    case 0:
        at.name = cell->text().toStdString();
        at.hidden = cell->checkState()/2;
        break;
    default:
        at.coord[column-1] = locale().toDouble(cell->text());
        at.fix[column-1] = cell->checkState()/2;
    }
    emit stepChanged();
}
