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
    ui->atomTableButton->setAttribute(Qt::WA_MacBrushedMetal, true);
    ui->cellWidgetButton->setAttribute(Qt::WA_MacBrushedMetal, true);
    ui->kpointStackButton->setAttribute(Qt::WA_MacBrushedMetal, true);
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
        ui->atomFmtBox->setCurrentIndex(static_cast<int>(curStep->getFmt()));
    }else if((change & Change::fmt) != 0){
        curStep = &master->curStep->asFmt(master->getFmt());
        ui->atomFmtBox->setCurrentIndex(static_cast<int>(master->getFmt()));
    }
    if ((change & (Change::atoms | Change::fmt)) != 0) {
        fillAtomTable();
    }
    if ((change & Change::cell) != 0) {
        fillCell();
    }
    if ((change & Change::kpoints) != 0) {
        fillKPoints();
    }
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
    if(!ui->atomTableButton->isChecked()){
        atomsOutdated = true;
        return;
    }
    //Fill atom list
    QSignalBlocker blockTable(ui->atomTable);
    int oldCount = ui->atomTable->rowCount();
    auto nat = static_cast<int>(curStep->getNat());
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
                    Qt::CheckState(static_cast<int>(at->properties[Hidden])*2));
        for(int k=0;k!=3;++k){
            ui->atomTable->item(j,k+1)->setText(QString::number(at->coord[k]));
            ui->atomTable->item(j,k+1)->setCheckState(
                        Qt::CheckState(at->properties[k]*2));
        }
        ++at;
    }
    atomsOutdated = false;
}

void MolWidget::fillKPoints()
{
    auto& kpoints = curMol->getKPoints();
    ui->activeKpoint->setCurrentIndex(static_cast<int>(kpoints.active));
    // fill mpg
    ui->mpg_x->setValue(kpoints.mpg.x);
    ui->mpg_y->setValue(kpoints.mpg.y);
    ui->mpg_z->setValue(kpoints.mpg.z);
    ui->mpg_x_off->setValue(static_cast<double>(kpoints.mpg.sx));
    ui->mpg_y_off->setValue(static_cast<double>(kpoints.mpg.sy));
    ui->mpg_z_off->setValue(static_cast<double>(kpoints.mpg.sz));
    // fill discrete
    auto discToCheckstate = [](const KPoints::Discrete&k, KPoints::Discrete::Properties p){
        if((k.properties&p) != 0){
            return Qt::CheckState::Checked;
        }
        return Qt::CheckState::Unchecked;
    };
    ui->crystal->setCheckState(discToCheckstate(kpoints.discrete, kpoints.discrete.crystal));
    ui->bands->setCheckState(discToCheckstate(kpoints.discrete, kpoints.discrete.band));
    auto& discretetable = *(ui->discretetable);
    int oldCount = discretetable.rowCount();
    auto newCount = static_cast<int>(kpoints.discrete.kpoints.size());
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
    vec[static_cast<size_t>(row)][static_cast<size_t>(column)] =
            locale().toFloat(ui->cellVecTable->item(row,column)->text());
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
    Atom at = (*curStep)[static_cast<size_t>(row)];
    const QTableWidgetItem *cell = ui->atomTable->item(row,column);
    bool checkState = (cell->checkState()/2) != 0;
    if (column == 0){
        if(at.properties[Hidden] != checkState){
            at.properties[Hidden] = checkState;
        }else{
            at.name = cell->text().toStdString();
        }
        triggerUpdate(Change::atoms);
    } else {
        const auto col = static_cast<size_t>(column-1);
        if(at.properties[col] != checkState){
            at.properties[col] = checkState;
        }else{
            at.coord[col] = locale().toFloat(cell->text());
            triggerUpdate(Change::atoms);
        }
    }
}

void MolWidget::on_atomFmtBox_currentIndexChanged(int index)
{
    master->setFmt(index, false, false);
}

void MolWidget::on_atomFmtButton_clicked()
{
    master->setFmt(ui->atomFmtBox->currentIndex(), true, false);
}

void MolWidget::on_molList_currentIndexChanged(int index)
{
    master->setMol(index);
}

void MolWidget::registerMol(const std::string& name)
{
    ui->molList->addItem(name.c_str());
    ui->molList->setCurrentIndex(ui->molList->count()-1);
}

void MolWidget::on_atomTableButton_toggled(bool checked)
{
    ui->atomContainer->setVisible(checked);
    if(checked && atomsOutdated){
        fillAtomTable();
    }
}
