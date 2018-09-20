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
    master->curStep->modWrap();
    triggerUpdate(GuiChange::atoms);
}

void CellModWidget::on_cropButton_clicked()
{
    master->curStep->modCrop();
    triggerUpdate(GuiChange::atoms);
}

void CellModWidget::on_multButton_clicked()
{
    master->curStep->modMultiply(
                static_cast<size_t>(ui->xMultSel->value()),
                static_cast<size_t>(ui->yMultSel->value()),
                static_cast<size_t>(ui->zMultSel->value()));
    triggerUpdate(GuiChange::atoms | GuiChange::cell);
}

void CellModWidget::on_alignButton_clicked()
{
    master->curStep->modAlign(
                static_cast<uint8_t>(ui->stepVecSel->currentIndex()),
                static_cast<uint8_t>(ui->coordVecSel->currentIndex()));
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
    master->curStep->modReshape(newMat,
                                ui->cdmSel->value(),
                                static_cast<CdmFmt>(ui->cdmFmtSel->currentIndex()));
    triggerUpdate(GuiChange::atoms | GuiChange::cell);
}
