#include "molwidget.h"
#include "ui_molwidget.h"
#include "molwidget_aux/newelement.h"
#include "vipsterapplication.h"
#include <QTableWidgetItem>
#include <QMessageBox>
#include <QMenu>

using namespace Vipster;


constexpr const char* inactiveKpoints[] = {"Gamma", "Monkhorst-Pack grid", "Discrete"};
constexpr const char* activeKpoints[] = {"Gamma (active)",
                                         "Monkhorst-Pack grid (active)",
                                         "Discrete (active)"};

MolWidget::MolWidget(QWidget *parent) :
    BaseWidget(parent),
    ui(new Ui::MolWidget)
{
    ui->setupUi(this);

    // Connect ui elements
    connect(ui->typeButton, &QPushButton::toggled, ui->typeContainer, &QWidget::setVisible);
    connect(ui->kpointButton, &QPushButton::toggled, ui->kpointContainer, &QWidget::setVisible);

    // TODO: hide all by default?
    ui->typeContainer->setVisible(false);

    // setup k-points
    ui->discretetable->addAction(ui->actionNew_K_Point);
    ui->discretetable->addAction(ui->actionDelete_K_Point);
    ui->kpointContainer->setVisible(ui->kpointButton->isChecked());

    // Connect to app state changes
    connect(&vApp, &Application::stepChanged, this, &MolWidget::updateStep);

    connect(&vApp, &Application::activeMolChanged, this, &MolWidget::setActiveMol);
}

MolWidget::~MolWidget()
{
    delete ui;
}

void MolWidget::updateStep(Step &step)
{
    // TODO: missing cell stuff?
}

void MolWidget::setActiveMol(Molecule &mol)
{
    //curMol = &mol;
    ui->activeKpoint->setCurrentIndex(static_cast<int>(mol.kpoints.active));
    fillKPoints();
}

void MolWidget::updateWidget(GUI::change_t change)
{
    if (!updateTriggered){
        if ((change & GUI::molChanged) == GUI::molChanged) {
      //      curMol = vApp.curMol;
        }
        if ((change & GUI::stepChanged) == GUI::stepChanged) {
//            setActiveStep(*vApp.curStep);
        }else if (change & (GUI::Change::atoms | GUI::Change::fmt)) {
//            atomModel.setStep(ownStep.get());
//            setSelection();
        }else if (change & (GUI::Change::selection)){
//            setSelection();
        }
    }
    if (change & GUI::Change::atoms) {
        ui->typeWidget->setTable(&vApp.curMol->getPTE());
        //checkOverlap();
        //bondModel.setStep(&*ownStep, vApp.stepdata[vApp.curStep].automatic_bonds);
    }
    if (change & GUI::Change::cell) {
        //fillCell();
    }
    if (change & GUI::Change::kpoints) {
        //ui->activeKpoint->setCurrentIndex(static_cast<int>(curMol->kpoints.active));
        //fillKPoints();
    }
}

void MolWidget::fillKPoints()
{
    QSignalBlocker blockCrystal(ui->crystal);
    QSignalBlocker blockBands(ui->bands);
    QSignalBlocker blockDisc(ui->discretetable);
    const auto& kpoints = vApp.curMol->kpoints;
    for(int i=0; i<3; ++i){
        if(i == static_cast<int>(kpoints.active)){
            ui->activeKpoint->setItemText(i, activeKpoints[i]);
        }else{
            ui->activeKpoint->setItemText(i, inactiveKpoints[i]);
        }
    }
    // fill mpg
    ui->mpg_x->setValue(kpoints.mpg.x);
    ui->mpg_y->setValue(kpoints.mpg.y);
    ui->mpg_z->setValue(kpoints.mpg.z);
    ui->mpg_x_off->setValue(static_cast<double>(kpoints.mpg.sx));
    ui->mpg_y_off->setValue(static_cast<double>(kpoints.mpg.sy));
    ui->mpg_z_off->setValue(static_cast<double>(kpoints.mpg.sz));
    // fill discrete
    ui->crystal->setCheckState((kpoints.discrete.properties & kpoints.discrete.crystal) ?
                                   Qt::CheckState::Checked :
                                   Qt::CheckState::Unchecked);
    ui->bands->setCheckState((kpoints.discrete.properties & kpoints.discrete.band) ?
                                   Qt::CheckState::Checked :
                                   Qt::CheckState::Unchecked);
    auto& discretetable = *(ui->discretetable);
    const auto& discpoints = kpoints.discrete.kpoints;
    auto count = static_cast<int>(discpoints.size());
    discretetable.clear();
    discretetable.setRowCount(count);
    for(int i=0; i<count; ++i){
        const auto& kp = discpoints[i];
        for(int j=0; j<3; ++j){
            discretetable.setItem(i,j, new QTableWidgetItem(QString::number(kp.pos[j])));
        }
        discretetable.setItem(i, 3, new QTableWidgetItem(QString::number(kp.weight)));
    }
}

void MolWidget::on_kFmtButton_clicked()
{
    auto oldFmt = static_cast<int>(vApp.curMol->kpoints.active);
    ui->activeKpoint->setItemText(oldFmt, inactiveKpoints[oldFmt]);
    auto newFmt = ui->activeKpoint->currentIndex();
    ui->activeKpoint->setItemText(newFmt, activeKpoints[newFmt]);
    vApp.curMol->kpoints.active = static_cast<KPoints::Fmt>(newFmt);
    triggerUpdate(GUI::Change::kpoints);
}

void MolWidget::on_bands_stateChanged(int arg)
{
    if(arg){
        vApp.curMol->kpoints.discrete.properties |= KPoints::Discrete::band;
    }else{
        vApp.curMol->kpoints.discrete.properties ^= KPoints::Discrete::band;
    }
    triggerUpdate(GUI::Change::kpoints);
}

void MolWidget::on_crystal_stateChanged(int arg)
{
    if(arg){
        vApp.curMol->kpoints.discrete.properties |= KPoints::Discrete::crystal;
    }else{
        vApp.curMol->kpoints.discrete.properties ^= KPoints::Discrete::crystal;
    }
    triggerUpdate(GUI::Change::kpoints);
}

void MolWidget::mpg_change()
{
    auto& kpoints = vApp.curMol->kpoints.mpg;
    if(sender() == ui->mpg_x){
        kpoints.x = ui->mpg_x->value();
    }else if(sender() == ui->mpg_y){
        kpoints.y = ui->mpg_y->value();
    }else if(sender() == ui->mpg_z){
        kpoints.z = ui->mpg_z->value();
    }else if(sender() == ui->mpg_x_off){
        kpoints.sx = ui->mpg_x_off->value();
    }else if(sender() == ui->mpg_y_off){
        kpoints.sy = ui->mpg_y_off->value();
    }else if(sender() == ui->mpg_z_off){
        kpoints.sz = ui->mpg_z_off->value();
    }
    triggerUpdate(GUI::Change::kpoints);
}

void MolWidget::on_discretetable_itemSelectionChanged()
{
    auto sel = ui->discretetable->selectedItems();
    if(sel.empty()){
        curKPoint = -1;
        ui->actionDelete_K_Point->setDisabled(true);
    }else{
        curKPoint = sel[0]->row();
        ui->actionDelete_K_Point->setEnabled(true);
    }
}

void MolWidget::on_actionNew_K_Point_triggered()
{
    auto& kpoints = vApp.curMol->kpoints.discrete.kpoints;
    kpoints.push_back(KPoints::Discrete::Point{});
    fillKPoints();
    triggerUpdate(GUI::Change::kpoints);
}

void MolWidget::on_actionDelete_K_Point_triggered()
{
    if(curKPoint < 0){
        throw Error{"MolWidget: \"Delete K-Point\" triggered with invalid selection"};
    }
    auto& kpoints = vApp.curMol->kpoints.discrete.kpoints;
    kpoints.erase(kpoints.begin()+curKPoint);
    fillKPoints();
    triggerUpdate(GUI::Change::kpoints);
}

void MolWidget::on_discretetable_cellChanged(int row, int column)
{
    auto& kp = vApp.curMol->kpoints.discrete.kpoints[row];
    QTableWidgetItem *cell = ui->discretetable->item(row, column);
    if(column == 3){
        kp.weight = cell->text().toDouble();
    }else{
        kp.pos[column] = cell->text().toDouble();
    }
    triggerUpdate(GUI::Change::kpoints);
}

void MolWidget::on_clearTableButton_clicked()
{
    vApp.curMol->cleanPTE();
    ui->typeWidget->setTable(&vApp.curMol->getPTE());
}

void MolWidget::on_newElemButton_clicked()
{
    if(newelement(vApp.curMol->getPTE()).exec()){
        ui->typeWidget->setTable(&vApp.curMol->getPTE());
    }
}
