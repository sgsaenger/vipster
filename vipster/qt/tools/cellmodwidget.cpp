#include "cellmodwidget.h"
#include "ui_cellmodwidget.h"

using namespace Vipster;

CellModWidget::CellModWidget(QWidget *parent) :
    BaseWidget(parent),
    ui(new Ui::CellModWidget)
{
    ui->setupUi(this);
}

CellModWidget::~CellModWidget()
{
    delete ui;
}

void CellModWidget::on_wrapButton_clicked()
{
    if(ui->trajecCheck->isChecked()){
        for(auto& step: master->curMol->getSteps()){
            step.modWrap();
        }
    }else{
        master->curStep->modWrap();
    }
    triggerUpdate(GuiChange::atoms);
}

void CellModWidget::on_cropButton_clicked()
{
    if(ui->trajecCheck->isChecked()){
        for(auto& step: master->curMol->getSteps()){
            step.modCrop();
        }
    }else{
        master->curStep->modCrop();
    }
    triggerUpdate(GuiChange::atoms);
}

void CellModWidget::on_multButton_clicked()
{
    auto x = static_cast<size_t>(ui->xMultSel->value());
    auto y = static_cast<size_t>(ui->yMultSel->value());
    auto z = static_cast<size_t>(ui->zMultSel->value());
    if(ui->trajecCheck->isChecked()){
        for(auto& step: master->curMol->getSteps()){
            step.modMultiply(x, y, z);
        }
    }else{
        master->curStep->modMultiply(x, y, z);
    }
    triggerUpdate(GuiChange::atoms | GuiChange::cell);
}

void CellModWidget::on_alignButton_clicked()
{
    auto stepdir = static_cast<uint8_t>(ui->stepVecSel->currentIndex());
    auto targetdir = static_cast<uint8_t>(ui->coordVecSel->currentIndex());
    if(ui->trajecCheck->isChecked()){
        for(auto& step: master->curMol->getSteps()){
            step.modAlign(stepdir, targetdir);
        }
    }else{
        master->curStep->modAlign(stepdir, targetdir);
    }
    triggerUpdate(GuiChange::atoms | GuiChange::cell);
}

void CellModWidget::on_reshapeButton_clicked()
{
    auto* table = ui->reshapeTable;
    Mat newMat;
    for(int i=0; i<3; ++i){
        for(int j=0; j<3; ++j){
            newMat[i][j] = table->item(i,j)->text().toFloat();
        }
    }
    auto cdm = static_cast<float>(ui->cdmSel->value());
    auto fmt = static_cast<CdmFmt>(ui->cdmFmtSel->currentIndex());
    if(ui->trajecCheck->isChecked()){
        for(auto& step: master->curMol->getSteps()){
            step.modReshape(newMat, cdm, fmt);
        }
    }else{
        master->curStep->modReshape(newMat, cdm, fmt);
    }
    triggerUpdate(GuiChange::atoms | GuiChange::cell);
}
