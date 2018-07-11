#include "data3dwidget.h"
#include "ui_data3dwidget.h"

using namespace Vipster;

Data3DWidget::Data3DWidget(QWidget *parent) :
    DataBase(parent),
    ui(new Ui::Data3DWidget)
{
    ui->setupUi(this);
}

Data3DWidget::~Data3DWidget()
{
    delete ui;
}

void Data3DWidget::setData(const BaseData* data)
{
    curData = dynamic_cast<const DataGrid3D_f*>(data);
    if(curData == nullptr){
        throw Error("Invalid dataset");
    }
    //init state
    auto dir = ui->sliceDir->currentIndex();
    ui->sliceVal->setValue(0);
    ui->sliceVal->setMaximum(curData->extent[dir]);
    ui->surfToggle->setCheckState(Qt::CheckState::Unchecked);
    //TODO: get min/max values from grid
//        ui->surfVal->setMinimum()
//        ui->surfVal->setMaximum()
}

void Data3DWidget::on_sliceBut_clicked()
{
    //TODO: toggle display
}

void Data3DWidget::on_sliceDir_currentIndexChanged(int index)
{
    ui->sliceVal->setValue(0);
    ui->sliceVal->setMaximum(curData->extent[index]);
}

void Data3DWidget::on_sliceVal_valueChanged(int)
{
    //TODO: create new plane, or toggle creation
}

void Data3DWidget::on_surfToggle_stateChanged(int arg1)
{

}

void Data3DWidget::on_surfVal_valueChanged(double arg1)
{

}

void Data3DWidget::on_surfBut_clicked()
{

}
