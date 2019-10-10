#include "data3dvwidget.h"
#include "ui_data3dvwidget.h"
#include "mainwindow.h"

using namespace Vipster;

Data3DVWidget::Data3DVWidget(QWidget *parent) :
    DataBase(parent),
    ui(new Ui::Data3DVWidget)
{
    ui->setupUi(this);
}

Data3DVWidget::~Data3DVWidget()
{
    delete ui;
}

void Data3DVWidget::setData(const BaseData* data)
{
    curData = dynamic_cast<const DataGrid3D_v*>(data);
    if(curData == nullptr){
        throw Error("Invalid dataset");
    }
}
