#include "../mainwindow.h"
#include "pinwidget.h"
#include "ui_pinwidget.h"
#include <QMessageBox>

using namespace Vipster;

PinWidget::PinWidget(QWidget *parent) :
    BaseWidget(parent),
    ui(new Ui::PinWidget)
{
    ui->setupUi(this);
}

PinWidget::~PinWidget()
{
    delete ui;
}

void PinWidget::updateWidget(GUI::change_t change)
{
//    if((change & GUI::vpChanged) == GUI::vpChanged){
//        // fill list with steps pinned in current viewport
//        QSignalBlocker block{this};
//        ui->stepList->clear();
//        curPin = nullptr;
//        pinnedSteps.clear();
//        for(auto &p : master->curVP->vpdata.extras){
//            auto test = std::dynamic_pointer_cast<PinnedStep>(p);
//            if(test){
//                pinnedSteps.push_back(std::move(test));
//                ui->stepList->addItem(pinnedSteps.back()->name.c_str());
//            }
//        }
    if((change & GUI::stepChanged) == GUI::stepChanged){
        // disable add-button when already pinned
        for(auto &dat: pinnedSteps){
            if(dat->curStep == master->curStep){
                ui->addStep->setDisabled(true);
                return;
            }
        }
        ui->addStep->setEnabled(true);
    }
    if(change & (GUI::Change::atoms | GUI::Change::trajec |
                       GUI::Change::cell | GUI::Change::settings)){
        // set gui-state
        on_stepList_currentRowChanged(ui->stepList->currentRow());
        // update GPU data
        for(auto &dat: pinnedSteps){
            const auto& settings = master->settings;
            dat->update(dat->curStep, settings.atRadVdW.val,
                        settings.atRadFac.val, settings.bondRad.val);
        }
    }
}

PinWidget::PinnedStep::PinnedStep(const GUI::GlobalData &glob,
                                  Step *step, const std::string& name,
                                  GUI::PBCVec mult)
    : GUI::StepData{glob, step},
      name{name},
      mult{mult}
{}

void PinWidget::PinnedStep::draw(const Vec &off,
                                 const GUI::PBCVec &m,
                                 const Mat &cv, bool drawCell, void *context)
{
    Vec off_loc = off + this->offset;
    Mat cv_loc = curStep->getCellVec() * curStep->getCellDim(CdmFmt::Bohr);
    auto mult_loc = curStep->hasCell() ? mult : GUI::PBCVec{{1,1,1}};
    if(repeat){
        for(int x=0; x<m[0]; ++x){
            for(int y=0; y<m[1]; ++y){
                for(int z=0; z<m[2]; ++z){
                    StepData::draw(off_loc + x*cv[0] + y*cv[1] + z*cv[2],
                            mult_loc, cv_loc, drawCell & this->cell, context);
                }
            }
        }
    }else{
        StepData::draw(off_loc, mult_loc, cv_loc, drawCell & this->cell, context);
    }
}

void PinWidget::setMult(int i)
{
    if (!curPin) return;
    auto &mult = curPin->mult;
    if(sender() == ui->xMultBox){
        mult[0] = static_cast<uint8_t>(i);
    }else if(sender() == ui->yMultBox){
        mult[1] = static_cast<uint8_t>(i);
    }else if(sender() == ui->zMultBox){
        mult[2] = static_cast<uint8_t>(i);
    }
    triggerUpdate(GUI::Change::extra);
}

void PinWidget::setOffset(double d)
{
    if (!curPin) return;
    auto &off = curPin->offset;
    if(sender() == ui->xOffset){
        off[0] = static_cast<float>(d);
    }else if(sender() == ui->yOffset){
        off[1] = static_cast<float>(d);
    }else if(sender() == ui->zOffset){
        off[2] = static_cast<float>(d);
    }
    triggerUpdate(GUI::Change::extra);
}

void PinWidget::on_showCell_toggled(bool checked)
{
    if (!curPin) return;
    curPin->cell = checked;
    triggerUpdate(GUI::Change::extra);
}

void PinWidget::on_repeatStep_toggled(bool checked)
{
    if (!curPin) return;
    curPin->repeat = checked;
    triggerUpdate(GUI::Change::extra);
}

void PinWidget::on_delStep_clicked()
{
    // remove from viewports
    for(auto& vp: master->viewports){
        auto& vpdata = vp->vpdata.extras;
        auto pos = std::find_if(vpdata.begin(), vpdata.end(),
                                [&](const auto& p){return curPin.get() == p.get();});
        if(pos != vpdata.end()){
            vpdata.erase(pos);
        }
    }
    // remove local infos
    ui->insertStep->setDisabled(true);
    if(curPin->curStep == master->curStep){
        ui->addStep->setEnabled(true);
    }
    auto pos2 = std::find_if(pinnedSteps.begin(), pinnedSteps.end(),
                            [&](const auto& p){return curPin.get() == p.get();});
    pinnedSteps.erase(pos2);
    delete ui->stepList->takeItem(ui->stepList->currentRow());
    triggerUpdate(GUI::Change::extra);
}

void PinWidget::on_addStep_clicked()
{
    ui->addStep->setDisabled(true);
    // add to list of steps
    pinnedSteps.push_back(std::make_shared<PinnedStep>(master->globals, master->curStep,
        master->curMol->getName() + " (Step " + std::to_string(master->curVP->moldata[master->curMol].curStep) + ')',
        GUI::PBCVec{1,1,1}));
    pinnedSteps.back()->update(pinnedSteps.back()->curStep,
                               master->settings.atRadVdW.val, master->settings.atRadFac.val,
                               master->settings.bondRad.val);
    ui->stepList->addItem(QString::fromStdString(pinnedSteps.back()->name));
    // enable in current viewport
    master->curVP->vpdata.extras.push_back(pinnedSteps.back());
    triggerUpdate(GUI::Change::extra);
}

void PinWidget::on_stepList_currentRowChanged(int currentRow)
{
    curPin = currentRow < 0 ? nullptr : pinnedSteps[currentRow];
    auto hasPin = static_cast<bool>(curPin);
    auto hasCell = hasPin ? curPin->curStep->hasCell() : false;
    ui->insertStep->setEnabled(hasPin ? curPin->curStep != master->curStep : false);
    ui->delStep->setEnabled(hasPin);
    ui->showStep->setEnabled(hasPin);
    ui->repeatStep->setEnabled(hasPin);
    ui->xOffset->setEnabled(hasPin);
    ui->yOffset->setEnabled(hasPin);
    ui->zOffset->setEnabled(hasPin);
    ui->showCell->setEnabled(hasCell);
    ui->xMultBox->setEnabled(hasCell);
    ui->yMultBox->setEnabled(hasCell);
    ui->zMultBox->setEnabled(hasCell);
    ui->xFit->setEnabled(hasCell);
    ui->yFit->setEnabled(hasCell);
    ui->zFit->setEnabled(hasCell);
    if(hasPin){
        QSignalBlocker block{ui->showStep};
        auto &e = master->curVP->vpdata.extras;
        auto pos = std::find(e.begin(), e.end(), curPin);
        ui->showStep->setChecked(pos != e.end());
        QSignalBlocker block1{ui->showCell};
        ui->showCell->setChecked(curPin->cell);
        QSignalBlocker block2{ui->repeatStep};
        ui->repeatStep->setChecked(curPin->repeat);
        QSignalBlocker blockx{ui->xMultBox};
        ui->xMultBox->setValue(curPin->mult[0]);
        QSignalBlocker blocky{ui->yMultBox};
        ui->yMultBox->setValue(curPin->mult[1]);
        QSignalBlocker blockz{ui->zMultBox};
        ui->zMultBox->setValue(curPin->mult[2]);
    }
}

void PinWidget::on_insertStep_clicked()
{
    if (!curPin || (curPin->curStep == master->curStep)) return;
    Step s = *curPin->curStep;
    s.modShift(s.formatVec(curPin->offset, AtomFmt::Bohr, s.getFmt()));
    std::array<bool,3> fit = {ui->xFit->isChecked(),
                              ui->yFit->isChecked(),
                              ui->zFit->isChecked()};
    if (s.hasCell() && (curPin->mult != GUI::PBCVec{{1,1,1}})) {
        const auto &m = curPin->mult;
        s.modMultiply(m[0], m[1], m[2]);
    }
    if (s.hasCell() && (fit != std::array<bool,3>{{false, false, false}})){
        auto fac = master->curStep->getCellDim(CdmFmt::Bohr) /
                s.getCellDim(CdmFmt::Bohr);
        auto cell = s.getCellVec();
        const auto& target = master->curStep->getCellVec();
        if (fit[0]) cell[0] = target[0] * fac;
        if (fit[1]) cell[1] = target[1] * fac;
        if (fit[2]) cell[2] = target[2] * fac;
        s.setCellVec(cell, true);
    }
    master->curStep->newAtoms(s);
    // immediately hide pinned step
    auto &e = master->curVP->vpdata.extras;
    auto pos = std::find(e.begin(), e.end(), curPin);
    if(pos != e.end()){
        e.erase(pos);
    }
    triggerUpdate(GUI::Change::atoms);
}

void PinWidget::on_showStep_toggled(bool checked)
{
    if (!curPin) return;
    if (checked) {
        // insert into viewports extras
        master->curVP->vpdata.extras.push_back(curPin);
    }else{
        // remove from viewports extras
        auto &e = master->curVP->vpdata.extras;
        auto pos = std::find(e.begin(), e.end(), curPin);
        if(pos != e.end()){
            e.erase(pos);
        }
    }
    triggerUpdate(GUI::Change::extra);
}

void PinWidget::on_helpButton_clicked()
{
    QMessageBox::information(this, QString("About pinning"),
        "Pinned Steps are drawn along the currently active Step.\n\n"
        "\"Repeating\" a Step means it is drawn in periodic repetitions "
        "of the active Step, i.e. with the periodicity of the active Step.\n"
        "Contrarily, specifying the multipliers for the pinned Step "
        "itself draws it with its own periodicity.\n"
        "Specifying the offset allows the pinned Step to be shifted against the active Step "
        "without having to modify its structure.\n\n"
        "Inserting a Step takes both offset and multipliers into account, "
        "and additionally performs cell vector fitting if requested.\n\n"
        "Cell vector fitting can be used e.g. to create commensurable super cells. "
        "If fitting is enabled for a direction, "
        "the lattice will be shrunken or stretched to "
        "match this dimension in the active Step."
    );
}
