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
    if ((change & molChanged) == molChanged) {
        curMol = master->curMol;
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
    ui->cellDimBox->setValue(static_cast<double>(
                                 curStep->getCellDim(
            static_cast<CdmFmt>(ui->cellFmt->currentIndex()))));
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
    auto& kpoints = curMol->getKPoints();
    ui->activeKpoint->setCurrentIndex(static_cast<int>(kpoints.active));
    // fill mpg
    ui->mpg_x->setValue(kpoints.mpg.x);
    ui->mpg_y->setValue(kpoints.mpg.y);
    ui->mpg_z->setValue(kpoints.mpg.z);
    ui->mpg_x_off->setValue(kpoints.mpg.sx);
    ui->mpg_y_off->setValue(kpoints.mpg.sy);
    ui->mpg_z_off->setValue(kpoints.mpg.sz);
    // fill discrete
    auto discToCheckstate = [](const KPoints::Discrete&k, KPoints::Discrete::Properties p){
        if(k.properties&p){
            return Qt::CheckState::Checked;
        }else{
            return Qt::CheckState::Unchecked;
        }
    };
    ui->crystal->setCheckState(discToCheckstate(kpoints.discrete, kpoints.discrete.crystal));
    ui->bands->setCheckState(discToCheckstate(kpoints.discrete, kpoints.discrete.band));
    auto& discretetable = *(ui->discretetable);
    int oldCount = discretetable.rowCount();
    int newCount = kpoints.discrete.kpoints.size();
    discretetable.setRowCount(newCount);
    if (oldCount < newCount) {
        for (int j=oldCount; j!=newCount; ++j) {
            for (int k=0; k!=4; ++k) {
                discretetable.setItem(j,k, new QTableWidgetItem());
            }
        }
    }
    auto kpoint = kpoints.discrete.kpoints.begin();
    for (int j=0; j!=newCount; ++j) {
        discretetable.item(j,0)->setText(QString::number(kpoint->pos[0]));
        discretetable.item(j,1)->setText(QString::number(kpoint->pos[1]));
        discretetable.item(j,2)->setText(QString::number(kpoint->pos[2]));
        discretetable.item(j,3)->setText(QString::number(kpoint->weight));
    }
}

void MolWidget::on_cellEnabled_toggled(bool checked)
{
    curStep->enableCell(checked);
    triggerUpdate(Change::cell);
}

void MolWidget::on_cellFmt_currentIndexChanged(int idx)
{
    QSignalBlocker blockCDB(ui->cellDimBox);
    ui->cellDimBox->setValue(static_cast<double>(curStep->getCellDim(static_cast<CdmFmt>(idx))));
}

void MolWidget::on_cellDimBox_valueChanged(double cdm)
{
    curStep->setCellDim(static_cast<float>(cdm),
                        static_cast<CdmFmt>(ui->cellFmt->currentIndex()),
                        ui->cellScaleBox->isChecked());
    // if needed, trigger atom update
    Change change = Change::cell;
    if(ui->cellScaleBox->isChecked() != (curStep->getFmt()>=AtomFmt::Crystal)){
        change = static_cast<Change>(Change::cell | Change::atoms);
        fillAtomTable();
    }
    ui->cellEnabled->setCheckState(Qt::CheckState::Checked);
    triggerUpdate(change);
}

void MolWidget::on_cellVecTable_cellChanged(int row, int column)
{
    Mat vec;
    vec = curStep->getCellVec();
    vec[row][column] = locale().toDouble(ui->cellVecTable->item(row,column)->text());
    curStep->setCellVec(vec, ui->cellScaleBox->isChecked());
    // if needed, trigger atom update
    Change change = Change::cell;
    if(ui->cellScaleBox->isChecked() != (curStep->getFmt()==AtomFmt::Crystal)){
        change = static_cast<Change>(Change::cell | Change::atoms);
        fillAtomTable();
    }
    ui->cellEnabled->setCheckState(Qt::CheckState::Checked);
    triggerUpdate(change);
}

void MolWidget::on_atomTable_cellChanged(int row, int column)
{
    Atom at = (*curStep)[row];
    const QTableWidgetItem *cell = ui->atomTable->item(row,column);
    if (column == 0){
        at.name = cell->text().toStdString();
        at.properties[Hidden] = cell->checkState()/2;
    } else {
        // TODO: property assignment toggles pse-reevaluation in evaluateCache!
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
