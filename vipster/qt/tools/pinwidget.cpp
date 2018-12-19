#include "pinwidget.h"
#include "ui_pinwidget.h"

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

void PinWidget::updateWidget(guiChange_t change)
{
    if((change & guiStepChanged) == guiStepChanged){
        // enable drawing of step previously selected in mainwindow
        if(mainStep){
            auto& dat = stepMap.at(mainStep);
            dat.gpu_data.update(mainStep, settings.showBonds.val,
                                settings.showCell.val & dat.cell & mainStep->hasCell());
            if(dat.display){
                master->addExtraData(&dat.gpu_data);
            }
            mainStep = nullptr;
        }
        // if current step would be drawn twice, deactivate
        auto pos = stepMap.find(master->curStep);
        if(pos != stepMap.end()){
            master->delExtraData(&pos->second.gpu_data);
            mainStep = master->curStep;
            ui->addStep->setDisabled(true);
        }else{
            ui->addStep->setEnabled(true);
        }
    }
}

void PinWidget::on_showStep_toggled(bool checked)
{
    auto& dat = stepMap.at(activeStep);
    dat.display = checked;
    stepMap.at(activeStep).display = checked;
    if(activeStep != mainStep){
        if(checked){
            master->addExtraData(&dat.gpu_data);
        }else{
            master->delExtraData(&dat.gpu_data);
        }
    }
    triggerUpdate(GuiChange::extra);
}

void PinWidget::on_showCell_toggled(bool checked)
{
    auto& dat = stepMap.at(activeStep);
    dat.cell = checked;
    dat.gpu_data.update(activeStep, settings.showBonds.val,
                        checked & settings.showCell.val & activeStep->hasCell());
    triggerUpdate(GuiChange::extra);
}

void PinWidget::on_delStep_clicked()
{
    auto pos = ui->stepList->currentRow();
    stepList.erase(stepList.begin()+pos);
    if(activeStep != mainStep){
        master->delExtraData(&stepMap.at(activeStep).gpu_data);
    }
    stepMap.erase(activeStep);
    delete ui->stepList->takeItem(pos);
}

void PinWidget::on_addStep_clicked()
{
    ui->addStep->setDisabled(true);
    mainStep = master->curStep;
    stepMap.emplace(mainStep, PinnedStep{
            true,
            settings.showCell.val,
            GUI::StepData{master->getGLGlobals(),
                          mainStep}});
    stepList.push_back(mainStep);
    ui->stepList->addItem(QString::fromStdString(master->curMol->getName() + " (Step " +
                          std::to_string(master->moldata[master->curMol].curStep) + ')'));
}

void PinWidget::on_stepList_currentRowChanged(int currentRow)
{
    if(currentRow < 0){
        activeStep = nullptr;
        ui->delStep->setDisabled(true);
        ui->showCell->setDisabled(true);
        ui->showStep->setDisabled(true);
        ui->insertStep->setDisabled(true);
    }else{
        activeStep = stepList[static_cast<size_t>(currentRow)];
        auto& dat = stepMap.at(activeStep);
        QSignalBlocker block{this};
        ui->delStep->setEnabled(true);
        ui->showCell->setEnabled(true);
        ui->showCell->setChecked(dat.cell);
        ui->showStep->setEnabled(true);
        ui->showStep->setChecked(dat.display);
        ui->insertStep->setEnabled(true);
    }
}

void PinWidget::on_insertStep_clicked()
{
    if(activeStep == master->curStep) return;
    master->curStep->newAtoms(activeStep->asFmt(master->curStep->getFmt()).getAtoms());
    triggerUpdate(GuiChange::atoms);
}
