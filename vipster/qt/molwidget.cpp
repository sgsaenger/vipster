#include "molwidget.h"
#include "ui_molwidget.h"
#include "mainwindow.h"
#include <QTableWidgetItem>
#include <QMessageBox>
#include <QMenu>

using namespace Vipster;

MolWidget::MolWidget(QWidget *parent) :
    BaseWidget(parent),
    ui(new Ui::MolWidget)
{
    ui->setupUi(this);
    ui->atomTableButton->setAttribute(Qt::WA_MacBrushedMetal, true);
    ui->cellWidgetButton->setAttribute(Qt::WA_MacBrushedMetal, true);
    ui->kpointStackButton->setAttribute(Qt::WA_MacBrushedMetal, true);
    ui->discretetable->addAction(ui->actionNew_K_Point);
    ui->discretetable->addAction(ui->actionDelete_K_Point);
    QSignalBlocker tableBlocker(ui->cellVecTable);
    for(int j=0;j!=3;++j){
        for(int k=0;k!=3;++k){
             ui->cellVecTable->setItem(j,k,new QTableWidgetItem());
        }
    }
    ui->atomTable->setModel(&molModel);
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
        molModel.setColumns(i);
    };
    for(auto& action: headerActions){
        connect(action, &QAction::toggled, this, changeColumns);
    }
    ui->atomTable->horizontalHeader()->setContextMenuPolicy(Qt::ActionsContextMenu);
    ui->atomTable->horizontalHeader()->addActions(headerActions);
}

MolWidget::~MolWidget()
{
    delete ui;
}

void MolWidget::updateWidget(guiChange_t change)
{
    if (updateTriggered) {
        updateTriggered = false;
        return;
    }
    if ((change & guiMolChanged) == guiMolChanged) {
        curMol = master->curMol;
    }
    if ((change & guiStepChanged) == guiStepChanged) {
        QSignalBlocker blockAtFmt(ui->atomFmtBox);
        // reset old fmt-string
        auto oldFmt = static_cast<int>(curStep.getFmt());
        ui->atomFmtBox->setItemText(oldFmt, inactiveFmt[oldFmt]);
        // assign StepFormatter to curStep, mark fmt as active
        auto fmt = master->curStep->getFmt();
        curStep = master->curStep->asFmt(fmt);
        molModel.setStep(&curStep);
        setSelection();
        auto ifmt = static_cast<int>(fmt);
        ui->atomFmtBox->setCurrentIndex(ifmt);
        ui->atomFmtBox->setItemText(ifmt, activeFmt[ifmt]);
    }else if (change & (GuiChange::atoms | GuiChange::fmt)) {
        molModel.setStep(&curStep);
        setSelection();
    }else if (change & (GuiChange::selection)){
        setSelection();
    }
    if (change & GuiChange::cell) {
        fillCell();
    }
    if (change & GuiChange::kpoints) {
        ui->activeKpoint->setCurrentIndex(static_cast<int>(curMol->getKPoints().active));
        fillKPoints();
    }
}

void MolWidget::fillCell()
{
    //Fill cell view
    QSignalBlocker blockCell(ui->cellVecTable);
    QSignalBlocker blockDim(ui->cellDimBox);
    QSignalBlocker blockEnabled(ui->cellEnabledBox);
    ui->cellEnabledBox->setChecked(curStep.hasCell());
    ui->cellDimBox->setValue(static_cast<double>(
                                 curStep.getCellDim(
            static_cast<CdmFmt>(ui->cellFmt->currentIndex()))));
    Mat vec = curStep.getCellVec();
    for(int j=0;j!=3;++j){
        for(int k=0;k!=3;++k){
            ui->cellVecTable->item(j,k)->setText(QString::number(vec[j][k]));
        }
    }
}

void MolWidget::on_cellTrajecButton_clicked()
{
    if(ui->cellEnabledBox->isChecked()){
        auto scale = ui->cellScaleBox->isChecked();
        auto dim = static_cast<float>(ui->cellDimBox->value());
        auto fmt = static_cast<CdmFmt>(ui->cellFmt->currentIndex());
        Mat vec{};
        for(int row=0; row<3; ++row){
            for(int col=0; col<3; ++col){
                vec[static_cast<size_t>(row)][static_cast<size_t>(col)] =
                    ui->cellVecTable->item(row,col)->text().toFloat();
            }
        }
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
    triggerUpdate(GuiChange::trajec);
}

void MolWidget::on_cellEnabledBox_toggled(bool checked)
{
    if(checked){
        guiChange_t change = GuiChange::cell;
        auto scale = ui->cellScaleBox->isChecked();
        if (scale) change |= GuiChange::atoms;
        auto dim = static_cast<float>(ui->cellDimBox->value());
        auto fmt = static_cast<CdmFmt>(ui->cellFmt->currentIndex());
        Mat vec{};
        for(int row=0; row<3; ++row){
            for(int col=0; col<3; ++col){
                vec[static_cast<size_t>(row)][static_cast<size_t>(col)] =
                    ui->cellVecTable->item(row,col)->text().toFloat();
            }
        }
        try{
            curStep.setCellVec(vec, scale);
        } catch(const Error& e){
            QMessageBox msg{this};
            msg.setText(QString{"Error setting cell vectors:\n"}+e.what());
            msg.exec();
            QSignalBlocker block{ui->cellEnabledBox};
            ui->cellEnabledBox->setCheckState(Qt::CheckState::Unchecked);
            return;
        }
        curStep.setCellDim(dim, fmt, scale);
        triggerUpdate(change);
    }else{
        curStep.enableCell(false);
        triggerUpdate(GuiChange::cell);
    }
}

void MolWidget::on_cellFmt_currentIndexChanged(int idx)
{
    QSignalBlocker blockCDB(ui->cellDimBox);
    ui->cellDimBox->setValue(static_cast<double>(curStep.getCellDim(static_cast<CdmFmt>(idx))));
}

void MolWidget::on_cellDimBox_valueChanged(double cdm)
{
    // if cell is disabled, exit early
    if(ui->cellEnabledBox->checkState() == Qt::CheckState::Unchecked){
        return;
    }
    auto dim = static_cast<float>(cdm);
    auto fmt = static_cast<CdmFmt>(ui->cellFmt->currentIndex());
    auto scale = ui->cellScaleBox->isChecked();
    curStep.setCellDim(dim, fmt, scale);
    guiChange_t change = GuiChange::cell;
    // if needed, trigger atom update
    if(scale){
        change |= GuiChange::atoms;
    }
    if(scale != (curStep.getFmt()>=AtomFmt::Crystal)){
        molModel.setStep(&curStep);
        setSelection();
    }
    triggerUpdate(change);
}

void MolWidget::on_cellVecTable_cellChanged(int row, int column)
{
    // if cell is disabled, exit early
    if(ui->cellEnabledBox->checkState() == Qt::CheckState::Unchecked){
        return;
    }
    Mat vec = curStep.getCellVec();
    vec[static_cast<size_t>(row)][static_cast<size_t>(column)] =
            ui->cellVecTable->item(row,column)->text().toFloat();
    auto scale = ui->cellScaleBox->isChecked();
    try{
        curStep.setCellVec(vec, scale);
    } catch(const Error& e){
        QMessageBox msg{this};
        msg.setText(QString{"Error setting cell vectors:\n"}+e.what());
        msg.exec();
        fillCell();
        return;
    }
    guiChange_t change = GuiChange::cell;
    // if needed, trigger atom update
    if(scale){
        change |= GuiChange::atoms;
    }
    if(scale != (curStep.getFmt()==AtomFmt::Crystal)){
        molModel.setStep(&curStep);
        setSelection();
    }
    triggerUpdate(change);
}

AtomFmt MolWidget::getAtomFmt()
{
    return static_cast<AtomFmt>(ui->atomFmtBox->currentIndex());
}

CdmFmt MolWidget::getCellFmt()
{
    return static_cast<CdmFmt>(ui->cellFmt->currentIndex());
}

void MolWidget::on_atomFmtBox_currentIndexChanged(int index)
{
    curStep = curStep.asFmt(static_cast<AtomFmt>(index));
    molModel.setStep(&curStep);
    setSelection();
}

void MolWidget::on_atomFmtButton_clicked()
{
    auto ifmt = ui->atomFmtBox->currentIndex();
    auto fmt = static_cast<AtomFmt>(ifmt);
    auto oldFmt = static_cast<int>(master->curStep->getFmt());
    ui->atomFmtBox->setItemText(oldFmt, inactiveFmt[oldFmt]);
    ui->atomFmtBox->setItemText(ifmt, activeFmt[ifmt]);
    master->curStep->setFmt(fmt);
    master->curSel->setFmt(fmt);
    if((fmt >= AtomFmt::Crystal) && !curStep.hasCell()){
        ui->cellEnabledBox->setChecked(true);
    }
    triggerUpdate(GuiChange::fmt);
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

void MolWidget::atomSelectionChanged(const QItemSelection &, const QItemSelection &)
{
    auto idx = ui->atomTable->selectionModel()->selectedRows();
    SelectionFilter filter{};
    filter.mode = SelectionFilter::Mode::Index;
    for(const auto& i: idx){
        filter.indices.emplace(static_cast<size_t>(i.row()), std::vector{SizeVec{}});
    }
    master->curSel->setFilter(filter);
    triggerUpdate(GuiChange::selection);
}

void MolWidget::setSelection()
{
    auto& table = ui->atomTable;
    auto selMod = table->selectionModel();
    // TODO: this is se problem, isn't it?
//    QSignalBlocker tableBlocker{selMod};
//    selMod->blockSignals(true);
    disconnect(ui->atomTable->selectionModel(), &QItemSelectionModel::selectionChanged,
            this, &MolWidget::atomSelectionChanged);
    table->clearSelection();
    table->setSelectionMode(QAbstractItemView::MultiSelection);
    for(const auto& i:master->curSel->getIndices()){
        table->selectRow(static_cast<int>(i.first));
//        selMod->select()
    }
//    selMod->blockSignals(false);
    connect(ui->atomTable->selectionModel(), &QItemSelectionModel::selectionChanged,
            this, &MolWidget::atomSelectionChanged);
//    table->setSelectionModel(table->selectionModel());
    update();
//    table->update();
    table->setSelectionMode(QAbstractItemView::ExtendedSelection);
//    emit selMod->
//    emit table->
}

void MolWidget::fillKPoints()
{
    const auto& kpoints = curMol->getKPoints();
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
    auto oldFmt = static_cast<int>(curMol->getKPoints().active);
    ui->activeKpoint->setItemText(oldFmt, inactiveKpoints[oldFmt]);
    auto newFmt = ui->activeKpoint->currentIndex();
    ui->activeKpoint->setItemText(newFmt, activeKpoints[newFmt]);
    curMol->getKPoints().active = static_cast<KPoints::Fmt>(newFmt);
    triggerUpdate(GuiChange::kpoints);
}

void MolWidget::on_bands_stateChanged(int arg)
{
    if(arg){
        curMol->getKPoints().discrete.properties |= KPoints::Discrete::band;
    }else{
        curMol->getKPoints().discrete.properties ^= KPoints::Discrete::band;
    }
    triggerUpdate(GuiChange::kpoints);
}

void MolWidget::on_crystal_stateChanged(int arg)
{
    if(arg){
        curMol->getKPoints().discrete.properties |= KPoints::Discrete::crystal;
    }else{
        curMol->getKPoints().discrete.properties ^= KPoints::Discrete::crystal;
    }
    triggerUpdate(GuiChange::kpoints);
}

void MolWidget::mpg_change()
{
    auto& kpoints = curMol->getKPoints().mpg;
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
    triggerUpdate(GuiChange::kpoints);
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
    auto& kpoints = curMol->getKPoints().discrete.kpoints;
    kpoints.push_back(KPoints::Discrete::Point{});
    fillKPoints();
    triggerUpdate(GuiChange::kpoints);
}

void MolWidget::on_actionDelete_K_Point_triggered()
{
    if(curKPoint < 0){
        throw Error{"MolWidget: \"Delete K-Point\" triggered with invalid selection"};
    }
    auto& kpoints = curMol->getKPoints().discrete.kpoints;
    kpoints.erase(kpoints.begin()+curKPoint);
    fillKPoints();
    triggerUpdate(GuiChange::kpoints);
}

void MolWidget::on_discretetable_cellChanged(int row, int column)
{
    auto& kp = curMol->getKPoints().discrete.kpoints[row];
    QTableWidgetItem *cell = ui->discretetable->item(row, column);
    if(column == 3){
        kp.weight = cell->text().toFloat();
    }else{
        kp.pos[column] = cell->text().toFloat();
    }
    triggerUpdate(GuiChange::kpoints);
}
