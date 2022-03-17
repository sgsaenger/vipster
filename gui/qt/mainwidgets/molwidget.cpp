#include "molwidget.h"
#include "ui_molwidget.h"
#include "../mainwindow.h"
#include "molwidget_aux/bonddelegate.h"
#include "molwidget_aux/doubledelegate.h"
#include "molwidget_aux/newelement.h"
#include <QTableWidgetItem>
#include <QMessageBox>
#include <QMenu>

using namespace Vipster;

MolWidget::MolWidget(QWidget *parent) :
    BaseWidget(parent),
    ui(new Ui::MolWidget)
{
    ui->setupUi(this);

    ui->typeContainer->setVisible(false);

    // setup k-points
    ui->discretetable->addAction(ui->actionNew_K_Point);
    ui->discretetable->addAction(ui->actionDelete_K_Point);
    ui->kpointContainer->setVisible(ui->kpointButton->isChecked());

    // setup cell-table
    ui->cellVecTable->setModel(&cellModel);
    ui->cellVecTable->setItemDelegate(new DoubleDelegate{});

    // setup atom table
    ui->atomTable->setModel(&atomModel);
    connect(ui->atomTable->selectionModel(), &QItemSelectionModel::selectionChanged,
            this, &MolWidget::atomSelectionChanged);
    headerActions.push_back(new QAction{"Show type", ui->atomTable});
    headerActions.push_back(new QAction{"Show coordinates", ui->atomTable});
    headerActions.push_back(new QAction{"Show charges", ui->atomTable});
    headerActions.push_back(new QAction{"Show forces", ui->atomTable});
    headerActions.push_back(new QAction{"Show visibility", ui->atomTable});
    headerActions.push_back(new QAction{"Show constraints", ui->atomTable});
    for(auto& action: headerActions){
        action->setCheckable(true);
    }
    headerActions[0]->setChecked(true);
    headerActions[1]->setChecked(true);
    auto changeColumns = [&](){
        // triggered through changed columns
        int i=0;
        for(int j=0; j<headerActions.size(); ++j){
            i += headerActions[j]->isChecked() << j;
        }
        atomModel.setColumns(i);
    };
    for(auto& action: headerActions){
        connect(action, &QAction::toggled, this, changeColumns);
    }
    ui->atomTable->horizontalHeader()->setContextMenuPolicy(Qt::ActionsContextMenu);
    ui->atomTable->horizontalHeader()->addActions(headerActions);
    ui->atomTable->setItemDelegate(new DoubleDelegate{});

    // setup bond table
    ui->bondTable->setModel(&bondModel);
    ui->bondTable->setItemDelegateForColumn(3, new BondDelegate{});
    ui->bondContainer->setVisible(ui->bondButton->isChecked());
}

MolWidget::~MolWidget()
{
    delete ui;
}

void MolWidget::updateWidget(GUI::change_t change)
{
    if (!updateTriggered){
        if ((change & GUI::molChanged) == GUI::molChanged) {
            curMol = master->curMol;
        }
        if ((change & GUI::stepChanged) == GUI::stepChanged) {
            // reset old fmt-string
            if(ownStep){
                auto oldFmt = static_cast<int>(ownStep->getFmt())+2;
                ui->atomFmtBox->setItemText(oldFmt, inactiveFmt[oldFmt]);
            }
            // assign StepFormatter to curStep, mark fmt as active
            auto fmt = master->curStep->getFmt();
            curStep = master->curStep;
            ownStep = std::make_unique<Step::formatter>(curStep->asFmt(fmt));
            atomModel.setStep(ownStep.get());
            setSelection();
            auto ifmt = static_cast<int>(fmt)+2;
            QSignalBlocker blockAtFmt(ui->atomFmtBox);
            ui->atomFmtBox->setCurrentIndex(ifmt);
            ui->atomFmtBox->setItemText(ifmt, activeFmt[ifmt]);
            // expose BondMode
            QSignalBlocker blockBondMode(ui->bondModeBox);
            bool autobonds = master->stepdata[curStep].automatic_bonds;
            ui->bondModeBox->setCurrentIndex(autobonds ? 1 : 0);
            if(autobonds){
                ui->bondSetButton->setDisabled(true);
            }else{
                ui->bondButton->setChecked(true);
                ui->bondSetButton->setEnabled(true);
            }
        }else if (change & (GUI::Change::atoms | GUI::Change::fmt)) {
            atomModel.setStep(ownStep.get());
            setSelection();
        }else if (change & (GUI::Change::selection)){
            setSelection();
        }
    }
    if (change & GUI::Change::atoms) {
        ui->typeWidget->setTable(&curMol->getPTE());
        checkOverlap();
        bondModel.setStep(&*ownStep, master->stepdata[curStep].automatic_bonds);
    }
    if (change & GUI::Change::cell) {
        fillCell();
    }
    if (change & GUI::Change::kpoints) {
        ui->activeKpoint->setCurrentIndex(static_cast<int>(curMol->kpoints.active));
        fillKPoints();
    }
}

void MolWidget::checkOverlap()
{
    const auto& ovlp = curStep->getOverlaps();
    if(ovlp.empty()){
        ui->ovlpTable->hide();
        ui->ovlpLabel->hide();
        ui->bondButton->setText("Bonds");
    }else{
        auto& table = *ui->ovlpTable;
        QSignalBlocker tableblock{table};
        table.setCurrentCell(-1, -1);
        table.setRowCount(ovlp.size());
        for(size_t i=0; i<ovlp.size(); ++i){
            table.setItem(i, 0, new QTableWidgetItem{QString::number(ovlp[i].at1)});
            table.setItem(i, 1, new QTableWidgetItem{QString::number(ovlp[i].at2)});
        }
        table.show();
        ui->ovlpLabel->show();
        ui->bondButton->setText("Bonds (!)");
    }
}

void MolWidget::fillCell()
{
    //Fill cell view
    QSignalBlocker blockDim(ui->cellDimBox);
    QSignalBlocker blockEnabled(ui->cellEnabledBox);
    ui->cellEnabledBox->setChecked(ownStep->hasCell());
    ui->cellDimBox->setValue(static_cast<double>(
                                 ownStep->getCellDim(
            static_cast<AtomFmt>(ui->cellFmt->currentIndex()))));
    cellModel.setStep(curStep);
}

void MolWidget::on_cellTrajecButton_clicked()
{
    if(ui->cellEnabledBox->isChecked()){
        auto scale = ui->cellScaleBox->isChecked();
        auto dim = ui->cellDimBox->value();
        auto fmt = static_cast<AtomFmt>(ui->cellFmt->currentIndex());
        Mat vec = cellModel.getVec();
        for(auto& step: master->curMol->getSteps()){
            if (&step == master->curStep) continue;
            step.setCellDim(dim, fmt, scale);
            step.setCellVec(vec, scale);
        }
    }else{
        for(auto& step: master->curMol->getSteps()){
            if (&step == master->curStep) continue;
            step.enableCell(false);
        }
    }
    triggerUpdate(GUI::Change::trajec);
}

void MolWidget::on_cellEnabledBox_toggled(bool checked)
{
    if(checked){
        GUI::change_t change = GUI::Change::cell;
        auto scale = ui->cellScaleBox->isChecked();
        if (scale) change |= GUI::Change::atoms;
        auto dim = ui->cellDimBox->value();
        auto fmt = static_cast<AtomFmt>(ui->cellFmt->currentIndex());
        try{
            auto vec = cellModel.getVec();
            if(vec == Mat{}){
                curStep->enableCell(true);
            }else{
                curStep->setCellVec(vec, scale);
            }
        } catch(const Error& e){
            QMessageBox::critical(this, "Error setting cell vectors", e.what());
            QSignalBlocker block{ui->cellEnabledBox};
            ui->cellEnabledBox->setCheckState(Qt::CheckState::Unchecked);
            return;
        }
        curStep->setCellDim(dim, fmt, scale);
        triggerUpdate(change);
    }else{
        curStep->enableCell(false);
        triggerUpdate(GUI::Change::cell);
    }
}

void MolWidget::on_cellFmt_currentIndexChanged(int idx)
{
    QSignalBlocker blockCDB(ui->cellDimBox);
    ui->cellDimBox->setValue(static_cast<double>(ownStep->getCellDim(static_cast<AtomFmt>(idx))));
}

void MolWidget::on_cellDimBox_valueChanged(double cdm)
{
    // if cell is disabled, exit early
    if(!ui->cellEnabledBox->isChecked()){
        return;
    }
    auto fmt = static_cast<AtomFmt>(ui->cellFmt->currentIndex());
    auto scale = ui->cellScaleBox->isChecked();
    curStep->setCellDim(cdm, fmt, scale);
    GUI::change_t change = GUI::Change::cell;
    // if needed, trigger atom update
    if(scale){
        change |= GUI::Change::atoms;
    }
    // short-circuit resetting the molModel
    if(scale == atomFmtAbsolute(ownStep->getFmt())){
        atomModel.setStep(ownStep.get());
        setSelection();
    }
    triggerUpdate(change);
}

bool MolWidget::scale()
{
    return ui->cellScaleBox->isChecked();
}

void MolWidget::on_atomFmtBox_currentIndexChanged(int index)
{
    ownStep = std::make_unique<Step::formatter>(curStep->asFmt(static_cast<AtomFmt>(index-2)));
    atomModel.setStep(ownStep.get());
    bondModel.setStep(ownStep.get(), master->stepdata[curStep].automatic_bonds);
    cellModel.setStep(curStep);
    setSelection();
}

void MolWidget::on_atomFmtButton_clicked()
{
    auto ifmt = ui->atomFmtBox->currentIndex();
    auto fmt = static_cast<AtomFmt>(ifmt-2);
    auto oldFmt = static_cast<int>(curStep->getFmt())+2;
    ui->atomFmtBox->setItemText(oldFmt, inactiveFmt[oldFmt]);
    ui->atomFmtBox->setItemText(ifmt, activeFmt[ifmt]);
    curStep->setFmt(fmt); // (possibly) invalidates dependent objects
    // reset formatter
    ownStep = std::make_unique<Step::formatter>(curStep->asFmt(fmt));
    // reset models
    atomModel.setStep(ownStep.get());
    bondModel.setStep(ownStep.get(), master->stepdata[curStep].automatic_bonds);
    cellModel.setStep(curStep);
    // reset selection
    SelectionFilter filter{};
    filter.mode = SelectionFilter::Mode::Index;
    filter.indices = master->curSel->getAtoms().indices;
    *master->curSel = curStep->select(filter);
    if(atomFmtRelative(fmt)){
        ui->cellEnabledBox->setChecked(true);
    }
    triggerUpdate(GUI::Change::fmt);
}

void MolWidget::on_atomHelpButton_clicked()
{
    QMessageBox::information(this, QString("About atoms"), Vipster::AtomsAbout);
}

void MolWidget::atomSelectionChanged(const QItemSelection &, const QItemSelection &)
{
    auto idx = ui->atomTable->selectionModel()->selectedRows();
    SelectionFilter filter{};
    filter.mode = SelectionFilter::Mode::Index;
    for(const auto& i: idx){
        filter.indices.emplace_back(static_cast<size_t>(i.row()), SizeVec{});
    }
    *master->curSel = curStep->select(filter);
    triggerUpdate(GUI::Change::selection);
}

void MolWidget::setSelection()
{
    auto& table = ui->atomTable;
    disconnect(ui->atomTable->selectionModel(), &QItemSelectionModel::selectionChanged,
            this, &MolWidget::atomSelectionChanged);
    table->clearSelection();
    table->setSelectionMode(QAbstractItemView::MultiSelection);
    for(const auto& i:master->curSel->getAtoms().indices){
        table->selectRow(static_cast<int>(i.first));
    }
    connect(ui->atomTable->selectionModel(), &QItemSelectionModel::selectionChanged,
            this, &MolWidget::atomSelectionChanged);
    update();
    table->setSelectionMode(QAbstractItemView::ExtendedSelection);
}

void MolWidget::fillKPoints()
{
    QSignalBlocker blockCrystal(ui->crystal);
    QSignalBlocker blockBands(ui->bands);
    QSignalBlocker blockDisc(ui->discretetable);
    const auto& kpoints = curMol->kpoints;
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
    auto oldFmt = static_cast<int>(curMol->kpoints.active);
    ui->activeKpoint->setItemText(oldFmt, inactiveKpoints[oldFmt]);
    auto newFmt = ui->activeKpoint->currentIndex();
    ui->activeKpoint->setItemText(newFmt, activeKpoints[newFmt]);
    curMol->kpoints.active = static_cast<KPoints::Fmt>(newFmt);
    triggerUpdate(GUI::Change::kpoints);
}

void MolWidget::on_bands_stateChanged(int arg)
{
    if(arg){
        curMol->kpoints.discrete.properties |= KPoints::Discrete::band;
    }else{
        curMol->kpoints.discrete.properties ^= KPoints::Discrete::band;
    }
    triggerUpdate(GUI::Change::kpoints);
}

void MolWidget::on_crystal_stateChanged(int arg)
{
    if(arg){
        curMol->kpoints.discrete.properties |= KPoints::Discrete::crystal;
    }else{
        curMol->kpoints.discrete.properties ^= KPoints::Discrete::crystal;
    }
    triggerUpdate(GUI::Change::kpoints);
}

void MolWidget::mpg_change()
{
    auto& kpoints = curMol->kpoints.mpg;
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
    auto& kpoints = curMol->kpoints.discrete.kpoints;
    kpoints.push_back(KPoints::Discrete::Point{});
    fillKPoints();
    triggerUpdate(GUI::Change::kpoints);
}

void MolWidget::on_actionDelete_K_Point_triggered()
{
    if(curKPoint < 0){
        throw Error{"MolWidget: \"Delete K-Point\" triggered with invalid selection"};
    }
    auto& kpoints = curMol->kpoints.discrete.kpoints;
    kpoints.erase(kpoints.begin()+curKPoint);
    fillKPoints();
    triggerUpdate(GUI::Change::kpoints);
}

void MolWidget::on_discretetable_cellChanged(int row, int column)
{
    auto& kp = curMol->kpoints.discrete.kpoints[row];
    QTableWidgetItem *cell = ui->discretetable->item(row, column);
    if(column == 3){
        kp.weight = cell->text().toDouble();
    }else{
        kp.pos[column] = cell->text().toDouble();
    }
    triggerUpdate(GUI::Change::kpoints);
}

void MolWidget::on_bondSetButton_clicked()
{
    ownStep->generateBonds();
    bondModel.setStep(ownStep.get(), master->stepdata[curStep].automatic_bonds);
    triggerUpdate(GUI::Change::atoms);
}

void MolWidget::on_bondHelpButton_clicked()
{
    QMessageBox::information(this, QString("About bonds"), Vipster::BondsAbout);
}

void MolWidget::on_bondModeBox_currentIndexChanged(int index)
{
    auto automatic = master->stepdata[curStep].automatic_bonds = static_cast<bool>(index);
    if(automatic){
        ui->bondSetButton->setDisabled(true);
    }else{
        ui->bondSetButton->setEnabled(true);
    }
    for(auto& vp: master->viewports){
        if(vp->curStep == curStep)
            vp->setBondMode(automatic);
    }
    bondModel.setStep(ownStep.get(), automatic);
    triggerUpdate(GUI::Change::atoms);
}

void MolWidget::on_ovlpTable_itemSelectionChanged()
{
    auto selection = ui->ovlpTable->selectedItems();
    if(!selection.empty()){
        const auto& sel = *selection[0];
        const auto& ovlp = curStep->getOverlaps()[sel.row()];
        auto idx = sel.column() == 0 ? ovlp.at1 : ovlp.at2;
        ui->atomTable->selectRow(idx);
    }
}

void MolWidget::on_clearTableButton_clicked()
{
    curMol->cleanPTE();
    ui->typeWidget->setTable(&curMol->getPTE());
}

void MolWidget::on_newElemButton_clicked()
{
    if(newelement(curMol->getPTE()).exec()){
        ui->typeWidget->setTable(&curMol->getPTE());
    }
}
