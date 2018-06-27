#include "cellmodwidget.h"
#include "ui_cellmodwidget.h"

using namespace Vipster;

CellModWidget::CellModWidget(QWidget *parent) :
    QWidget(parent),
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
    triggerUpdate(Change::atoms);
}

void CellModWidget::on_cropButton_clicked()
{
    master->curStep->modCrop();
    triggerUpdate(Change::atoms);
}

void CellModWidget::on_multButton_clicked()
{
    master->curStep->modMultiply(
                ui->xMultSel->value(),
                ui->yMultSel->value(),
                ui->zMultSel->value());
    triggerUpdate(Change::atoms);
}

void CellModWidget::on_alignButton_clicked()
{

}

void CellModWidget::on_reshapeButton_clicked()
{

}
