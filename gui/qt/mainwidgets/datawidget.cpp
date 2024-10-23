#include "datawidget.h"
#include "ui_datawidget.h"
#include "../mainwindow.h"

using namespace Vipster;

DataBase::DataBase(QWidget *parent)
    :QWidget{parent}
{}

DataWidget::DataWidget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::DataWidget)
{
    ui->setupUi(this);

    connect(&vApp, &MainWindow::dataListChanged, [&](){
        if (vApp.data.size() > ui->DataSel->count()) {
           const BaseData& dat = *vApp.data.back();
           ui->DataSel->addItem(dat.name.c_str());
           ui->DataSel->setCurrentIndex(ui->DataSel->count()-1);
        }
    });

    connect(ui->DataSel, &QComboBox::currentIndexChanged, this, &DataWidget::selectData);
}

DataWidget::~DataWidget()
{
    delete ui;
}

void DataWidget::selectData(int index)
{
    curData = std::next(vApp.data.begin(), index)->get();
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
