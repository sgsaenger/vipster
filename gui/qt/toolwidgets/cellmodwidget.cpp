#include "../mainwindow.h"
#include "cellmodwidget.h"
#include "ui_cellmodwidget.h"
#include "../vipsterapplication.h"
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

void CellModWidget::updateWidget(Vipster::GUI::change_t change)
{
    if(change & GUI::Change::cell){
        setEnabled(vApp.curStep().hasCell());
    }
}

void CellModWidget::on_wrapButton_clicked()
{
    if(ui->trajecCheck->isChecked()){
        vApp.invokeOnTrajec(&Step::modWrap);
    }else{
        vApp.invokeOnStep(&Step::modWrap);
    }
}

void CellModWidget::on_cropButton_clicked()
{
    if(ui->trajecCheck->isChecked()){
        vApp.invokeOnTrajec(&Step::modCrop);
    }else{
        vApp.invokeOnStep(&Step::modCrop);
    }
}

void CellModWidget::on_multButton_clicked()
{
    const auto x = static_cast<size_t>(ui->xMultSel->value());
    const auto y = static_cast<size_t>(ui->yMultSel->value());
    const auto z = static_cast<size_t>(ui->zMultSel->value());
    if(ui->trajecCheck->isChecked()){
        vApp.invokeOnTrajec(&Step::modMultiply, x, y, z);
    }else{
        vApp.invokeOnStep(&Step::modMultiply, x, y, z);
    }
}

void CellModWidget::on_alignButton_clicked()
{
    const auto stepdir = static_cast<uint8_t>(ui->stepVecSel->currentIndex());
    const auto targetdir = static_cast<uint8_t>(ui->coordVecSel->currentIndex());
    if(ui->trajecCheck->isChecked()){
        vApp.invokeOnTrajec(&Step::modAlign, stepdir, targetdir);
    }else{
        vApp.invokeOnStep(&Step::modAlign, stepdir, targetdir);
    }
}

void CellModWidget::on_reshapeButton_clicked()
{
    const auto* table = ui->reshapeTable;
    Mat newMat;
    for(int i=0; i<3; ++i){
        for(int j=0; j<3; ++j){
            newMat[i][j] = table->item(i,j)->text().toDouble();
        }
    }
    try{
        Mat_inv(newMat);
    }catch(Error& e){
        QMessageBox::critical(this, "Could not reshape cell", e.what());
        return;
    }
    const auto cdm = ui->cdmSel->value();
    const auto fmt = static_cast<AtomFmt>(ui->cdmFmtSel->currentIndex());
    if(ui->trajecCheck->isChecked()){
        vApp.invokeOnTrajec(&Step::modReshape, newMat, cdm, fmt);
    }else{
        vApp.invokeOnStep(&Step::modReshape, newMat, cdm, fmt);
    }
}
