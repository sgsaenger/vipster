#include "datawidget.h"
#include "ui_datawidget.h"
#include "mainwindow.h"

using namespace Vipster;

DataBase::DataBase(QWidget *parent)
    :BaseWidget{parent}
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

void DataWidget::updateWidget(GUI::change_t change)
{
    if(change & GUI::Change::data){
        if(master->data.size() > ui->DataSel->count()){
            const BaseData& dat = *master->data.back();
            ui->DataSel->addItem(dat.name.c_str());
            ui->DataSel->setCurrentIndex(ui->DataSel->count()-1);
        }
    }
    ui->ThreeDWidget->updateWidget(change);
    ui->VecWidget->updateWidget(change);
    ui->TwoDWidget->updateWidget(change);
}

void DataWidget::on_DataSel_currentIndexChanged(int index)
{
    curData = std::next(master->data.begin(), index)->get();
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
