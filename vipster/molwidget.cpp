#include "molwidget.h"
#include "ui_molwidget.h"
#include "atom.h"
#include <QTableWidgetItem>

using namespace Vipster;

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

void MolWidget::updateWidget(Change change)
{
    if (updateTriggered) {
        updateTriggered = false;
        return;
    }
    QSignalBlocker blockAtFmt(ui->atomFmtBox);
    if ((change & stepChanged) == stepChanged) {
        curStep = master->curStep;
        ui->atomFmtBox->setCurrentIndex((int)curStep->getFmt());
    }else if(change & Change::fmt){
        curStep = &master->curStep->asFmt(master->getFmt());
        ui->atomFmtBox->setCurrentIndex((int)master->getFmt());
    }
    if (change & (Change::atoms | Change::fmt))
        fillAtomTable();
    if (change & Change::cell)
        fillCell();
    if (change & Change::kpoints)
        fillKPoints();
}

void MolWidget::fillCell()
{
    //Fill cell view
    QSignalBlocker blockCell(ui->cellVecTable);
    QSignalBlocker blockDim(ui->cellDimBox);
    QSignalBlocker blockEnabled(ui->cellEnabled);
    ui->cellEnabled->setChecked(curStep->hasCell());
    ui->cellDimBox->setValue( curStep->getCellDim(
            (CdmFmt)ui->cellFmt->currentIndex()));
    Mat vec = curStep->getCellVec();
    for(int j=0;j!=3;++j){
        for(int k=0;k!=3;++k){
            ui->cellVecTable->item(j,k)->setText(QString::number(vec[j][k]));
        }
    }
}

void MolWidget::fillAtomTable(void)
{
    //Fill atom list
    QSignalBlocker blockTable(ui->atomTable);
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
    auto at = curStep->begin();
    for(int j=0;j!=nat;++j){
        //TODO: Fmt
        ui->atomTable->item(j,0)->setText(at->name.c_str());
        ui->atomTable->item(j,0)->setCheckState(
                    Qt::CheckState(at->properties[Hidden]*2));
        for(int k=0;k!=3;++k){
            ui->atomTable->item(j,k+1)->setText(QString::number(at->coord[k]));
            ui->atomTable->item(j,k+1)->setCheckState(
                        Qt::CheckState(at->properties[k]*2));
        }
        ++at;
    }
}

void MolWidget::fillKPoints()
{

}

void MolWidget::on_cellEnabled_toggled(bool checked)
{
    curStep->enableCell(checked);
    triggerUpdate(Change::cell);
}

void MolWidget::on_cellFmt_currentIndexChanged(int idx)
{
    QSignalBlocker blockCDB(ui->cellDimBox);
    ui->cellDimBox->setValue(curStep->getCellDim((CdmFmt)idx));
}

void MolWidget::on_cellDimBox_valueChanged(double cdm)
{
    curStep->setCellDim(cdm, (CdmFmt)ui->cellFmt->currentIndex(), ui->cellScaleBox->isChecked());
    triggerUpdate(Change::cell);
}

void MolWidget::on_cellVecTable_cellChanged(int row, int column)
{
    Mat vec;
    vec = curStep->getCellVec();
    vec[row][column] = locale().toDouble(ui->cellVecTable->item(row,column)->text());
    curStep->setCellVec(vec, ui->cellScaleBox->isChecked());
    triggerUpdate(Change::cell);
}

void MolWidget::on_atomTable_cellChanged(int row, int column)
{
    Atom at = (*curStep)[row];
    const QTableWidgetItem *cell = ui->atomTable->item(row,column);
    if (column == 0){
        at.name = cell->text().toStdString();
        at.properties[Hidden] = cell->checkState()/2;
    } else {
        at.coord[column-1] = locale().toDouble(cell->text());
        at.properties[column-1] = cell->checkState()/2;
    }
    triggerUpdate(Change::atoms);
}

void MolWidget::on_atomFmtBox_currentIndexChanged(int index)
{
    master->setFmt(index, false, false);
}

void MolWidget::on_atomFmtButton_clicked()
{
    master->setFmt(ui->atomFmtBox->currentIndex(), true, false);
}
