#include "datawidget.h"
#include "ui_datawidget.h"

using namespace Vipster;

DataBase::DataBase(QWidget *parent)
    :QWidget{parent}
{}

DataWidget::DataWidget(QWidget *parent) :
    BaseWidget(parent),
    ui(new Ui::DataWidget)
{
    ui->setupUi(this);
}

DataWidget::~DataWidget()
{
    delete ui;
}

void DataWidget::registerData(const std::string& name)
{
    ui->DataSel->addItem(name.c_str());
    ui->DataSel->setCurrentIndex(ui->DataSel->count()-1);
}

void DataWidget::on_DataSel_currentIndexChanged(int index)
{
    curData = master->data.at(index).get();
    if(dynamic_cast<const DataGrid3D_f*>(curData) != nullptr){
        ui->DataStack->setCurrentWidget(ui->ThreeDWidget);
        ui->ThreeDWidget->setData(curData);
    }else if(dynamic_cast<const DataGrid3D_v*>(curData) != nullptr){
        ui->DataStack->setCurrentWidget(ui->VecWidget);
        ui->VecWidget->setData(curData);
    }else if(dynamic_cast<const DataGrid2D_f*>(curData) != nullptr){
        ui->DataStack->setCurrentWidget(ui->TwoDWidget);
        ui->TwoDWidget->setData(curData);
    }else{
        throw Error("Invalid data set");
    }
}