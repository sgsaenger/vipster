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
    QSignalBlocker block{this};
    auto pos = planes.find(curData);
    if(pos != planes.end()){
        curPlane = &pos->second;
        ui->sliceBut->setChecked(curPlane->display);
    }else{
        curPlane = nullptr;
        ui->sliceBut->setChecked(false);
    }
}

void Data2DWidget::on_sliceBut_toggled(bool checked)
{
    if(curPlane){
        curPlane->display = checked;
        if(checked){
            master->addExtraData(&curPlane->gpu_data);
        }else{
            master->delExtraData(&curPlane->gpu_data);
        }
    }else if(checked){
        auto tmp = planes.emplace(curData, DatPlane{
            true,
            GUI::MeshData{master->getGLGlobals(),
                          {{{{0,0,0}},{{0,1,0}},{{1,1,0}},
                            {{0,0,0}},{{1,0,0}},{{1,1,0}}}},
                          curData->origin,
                          curData->cell,
                          settings.milCol.val}
            });
        curPlane = &tmp.first->second;
        master->addExtraData(&curPlane->gpu_data);
    }
}
