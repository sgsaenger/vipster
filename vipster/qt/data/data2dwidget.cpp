#include "data2dwidget.h"
#include "ui_data2dwidget.h"

using namespace Vipster;

Data2DWidget::Data2DWidget(QWidget *parent) :
    DataBase(parent),
    ui(new Ui::Data2DWidget)
{
    ui->setupUi(this);
}

Data2DWidget::~Data2DWidget()
{
    delete ui;
}

void Data2DWidget::setData(const BaseData* data)
{
    curData = dynamic_cast<const DataGrid2D_f*>(data);
    if(!curData){
        throw Error("Invalid dataset");
    }
    //TODO: toggle display?
    auto dir = ui->sliceDir->currentIndex();
    ui->sliceVal->setValue(0);
    ui->sliceVal->setMaximum(curData->extent[dir]);
}

void Data2DWidget::on_sliceBut_clicked()
{
    //TODO: toggle display
}

void Data2DWidget::on_sliceDir_currentIndexChanged(int index)
{
    ui->sliceVal->setValue(0);
    ui->sliceVal->setMaximum(curData->extent[index]);
}

void Data2DWidget::on_sliceVal_valueChanged(int)
{
    //TODO: create new plane, or toggle creation
}
