#include "data2dwidget.h"
#include "ui_data2dwidget.h"
#include "mainwindow.h"

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
        GUI::MeshData::Texture texture;
        texture.width = static_cast<int>(curData->extent[0]);
        texture.height = static_cast<int>(curData->extent[1]);
        auto minmax = std::minmax_element(curData->begin(), curData->end());
        auto min = *minmax.first;
        auto max = *minmax.second;
        auto factor = 100/(max-min);
        for(size_t y=0; y<curData->extent[1]; ++y){
            for(size_t x=0; x<curData->extent[0]; ++x){
                const auto& val = curData->operator()(x,y);
                auto tmp = (val-min)*factor;
                texture.data.push_back({static_cast<uint8_t>(std::round(2.55f*tmp)),
                                        static_cast<uint8_t>(std::round(2.55f*(100-abs(2*tmp-100)))),
                                        static_cast<uint8_t>(std::round(2.55f*(100-tmp))),
                                        128});
            }
        }
        auto tmp = planes.emplace(curData, DatPlane{
            true,
            GUI::MeshData{master->getGLGlobals(),
                          {{{0,0,0},{},{0,0}},{{0,1,0},{},{0,1}},{{1,1,0},{},{1,1}},
                           {{0,0,0},{},{0,0}},{{1,0,0},{},{1,0}},{{1,1,0},{},{1,1}}},
                          curData->origin,
                          curData->cell,
                          texture}
            });
        curPlane = &tmp.first->second;
        master->addExtraData(&curPlane->gpu_data);
    }
}
