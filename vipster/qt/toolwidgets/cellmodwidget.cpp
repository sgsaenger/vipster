#include "../mainwindow.h"
#include "cellmodwidget.h"
#include "ui_cellmodwidget.h"
#include <QMessageBox>

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
    GUI::change_t change;
    if(ui->trajecCheck->isChecked()){
        for(auto& step: master->curMol->getSteps()){
            step.modWrap();
        }
        change = GUI::Change::atoms | GUI::Change::trajec;
    }else{
        master->curStep->modWrap();
        change = GUI::Change::atoms;
    }
    triggerUpdate(change);
}

void CellModWidget::on_cropButton_clicked()
{
    GUI::change_t change;
    if(ui->trajecCheck->isChecked()){
        for(auto& step: master->curMol->getSteps()){
            step.modCrop();
        }
        change = GUI::Change::atoms | GUI::Change::trajec;
    }else{
        master->curStep->modCrop();
        change = GUI::Change::atoms;
    }
    triggerUpdate(change);
}

void CellModWidget::on_multButton_clicked()
{
    GUI::change_t change;
    auto x = static_cast<size_t>(ui->xMultSel->value());
    auto y = static_cast<size_t>(ui->yMultSel->value());
    auto z = static_cast<size_t>(ui->zMultSel->value());
    if(ui->trajecCheck->isChecked()){
        for(auto& step: master->curMol->getSteps()){
            step.modMultiply(x, y, z);
        }
        change = GUI::Change::atoms | GUI::Change::cell | GUI::Change::trajec;
    }else{
        master->curStep->modMultiply(x, y, z);
        change = GUI::Change::atoms | GUI::Change::cell;
    }
    triggerUpdate(change);
}

void CellModWidget::on_alignButton_clicked()
{
    auto stepdir = static_cast<uint8_t>(ui->stepVecSel->currentIndex());
    auto targetdir = static_cast<uint8_t>(ui->coordVecSel->currentIndex());
    GUI::change_t change;
    if(ui->trajecCheck->isChecked()){
        for(auto& step: master->curMol->getSteps()){
            step.modAlign(stepdir, targetdir);
        }
        change = GUI::Change::atoms | GUI::Change::cell | GUI::Change::trajec;
    }else{
        master->curStep->modAlign(stepdir, targetdir);
        change = GUI::Change::atoms | GUI::Change::cell;
    }
    triggerUpdate(change);
}

void CellModWidget::on_reshapeButton_clicked()
{
    GUI::change_t change;
    auto* table = ui->reshapeTable;
    Mat newMat;
    for(int i=0; i<3; ++i){
        for(int j=0; j<3; ++j){
            newMat[i][j] = table->item(i,j)->text().toFloat();
        }
    }
    try{
        Mat_inv(newMat);
    }catch(Error& e){
        QMessageBox::critical(this, "Could not reshape cell", e.what());
        return;
    }
    auto cdm = static_cast<float>(ui->cdmSel->value());
    auto fmt = static_cast<CdmFmt>(ui->cdmFmtSel->currentIndex());
    if(ui->trajecCheck->isChecked()){
        for(auto& step: master->curMol->getSteps()){
            step.modReshape(newMat, cdm, fmt);
        }
        change = GUI::Change::atoms | GUI::Change::cell | GUI::Change::trajec;
    }else{
        master->curStep->modReshape(newMat, cdm, fmt);
        change = GUI::Change::atoms | GUI::Change::cell;
    }
    triggerUpdate(change);
}
