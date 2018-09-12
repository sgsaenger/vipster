#include "data3dwidget.h"
#include "ui_data3dwidget.h"

using namespace Vipster;

Data3DWidget::Data3DWidget(QWidget *parent) :
    DataBase(parent),
    ui(new Ui::Data3DWidget)
{
    ui->setupUi(this);
    ui->surfVal->setValidator(&validator);
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
    QSignalBlocker block{this};

    // init plane-state
    auto slicePos = slices.find(curData);
    if(slicePos != slices.end()){
        curSlice = &slicePos->second;
        ui->sliceBut->setChecked(curSlice->display);
        ui->sliceDir->setCurrentIndex(static_cast<int>(curSlice->dir));
        ui->sliceVal->setMaximum(static_cast<int>(curData->extent[curSlice->dir]));
        ui->sliceVal->setValue(static_cast<int>(curSlice->pos));
    }else{
        curSlice = nullptr;
        ui->sliceBut->setChecked(false);
        ui->sliceDir->setCurrentIndex(0);
        ui->sliceVal->setMaximum(static_cast<int>(curData->extent[0]));
        ui->sliceVal->setValue(0);
    }

    // init surf-state
    auto minmax = std::minmax_element(curData->begin(), curData->end());
    auto min = static_cast<double>(*minmax.first);
    auto max = static_cast<double>(*minmax.second);
    auto precision = std::abs(max - min)/1000.;
    auto dec = -static_cast<int>(std::round(std::log10(precision)))+2;
    validator.setRange(min, max, dec);
    ui->maxLabel->setText(QString::number(max, 'g', dec));
    ui->minLabel->setText(QString::number(min, 'g', dec));
    auto surfPos = surfaces.find(curData);
    if(surfPos != surfaces.end()){
        curSurf = &surfPos->second;
        ui->surfBut->setChecked(curSurf->display);
        ui->surfToggle->setCheckState(Qt::CheckState(curSurf->plusmin*2));
        ui->surfVal->setText(QString::number(curSurf->isoval));
        auto _val = (curSurf->isoval - validator.bottom()) *
                    ui->surfSlider->maximum() /
                    (validator.top()-validator.bottom());
        ui->surfSlider->setValue(static_cast<int>(_val));
    }else{
        curSurf = nullptr;
        ui->surfBut->setChecked(false);
        ui->surfToggle->setCheckState(Qt::CheckState::Unchecked);
        ui->surfVal->setText("0.0");
        ui->surfSlider->setValue(0);
    }
}

std::vector<Vec> mkSlice(size_t dir)
{
    switch(dir){
    case 0:
        return {{{{0,0,0}},{{0,1,0}},{{0,1,1}},
                 {{0,0,0}},{{0,0,1}},{{0,1,1}}}};
    case 1:
        return {{{{0,0,0}},{{1,0,0}},{{1,0,1}},
                 {{0,0,0}},{{0,0,1}},{{1,0,1}}}};
    case 2:
        return {{{{0,0,0}},{{0,1,0}},{{1,1,0}},
                 {{0,0,0}},{{1,0,0}},{{1,1,0}}}};
    default:
        throw Error("Invalid direction for data slicing");
    }
}

void Data3DWidget::on_sliceDir_currentIndexChanged(int index)
{
    auto _index = static_cast<size_t>(index);
    //block sliceVal from triggering when sliceVal is higher than new max
    QSignalBlocker block{ui->sliceVal};
    ui->sliceVal->setMaximum(static_cast<int>(curData->extent[_index]));
    if(curSlice){
        curSlice->dir = _index;
        curSlice->gpu_data.update(mkSlice(_index));
        auto newPos = static_cast<int>(std::min(curSlice->pos, curData->extent[_index]));
        ui->sliceVal->setValue(newPos);
        //trigger manually
        on_sliceVal_valueChanged(newPos);
    }
}

void Data3DWidget::on_sliceVal_valueChanged(int pos)
{
    if(curSlice){
        curSlice->pos = static_cast<size_t>(pos);
        auto off = static_cast<float>(pos)/curData->extent[curSlice->dir];
        curSlice->gpu_data.update(curData->origin + off*curData->cell[curSlice->dir]);
        triggerUpdate(GuiChange::extra);
    }
}

void Data3DWidget::on_sliceBut_toggled(bool checked)
{
    if(curSlice){
        curSlice->display = checked;
        if(checked){
            master->addExtraData(&curSlice->gpu_data);
        }else{
            master->delExtraData(&curSlice->gpu_data);
        }
    }else if(checked){
        auto dir = static_cast<size_t>(ui->sliceDir->currentIndex());
        auto pos = static_cast<size_t>(ui->sliceVal->value());
        auto off = static_cast<float>(pos) / curData->extent[dir];
        auto tmp = slices.emplace(curData, DatSlice{
            true, dir, pos,
            GUI::MeshData{master->getGLGlobals(),
                          mkSlice(dir),
                          curData->origin + off*curData->cell[dir],
                          curData->cell,
                          settings.milCol.val}
            });
        curSlice = &tmp.first->second;
        master->addExtraData(&curSlice->gpu_data);
    }
}

void Data3DWidget::on_surfToggle_stateChanged(int state)
{
    if(curSurf){
        curSurf->plusmin = state;
    }
}

void Data3DWidget::on_surfSlider_valueChanged(int val)
{
    QSignalBlocker block{ui->surfVal};
    auto _val = (val * (validator.top()-validator.bottom()) /
                 ui->surfSlider->maximum()) + validator.bottom();
    ui->surfVal->setText(QString::number(_val));
    if(curSurf){
        curSurf->isoval = _val;
    }
}

void Data3DWidget::on_surfVal_editingFinished()
{
    QSignalBlocker block{ui->surfSlider};
    auto val = ui->surfVal->text().toDouble();
    auto _val = (val - validator.bottom()) *
                ui->surfSlider->maximum() /
                (validator.top()-validator.bottom());
    ui->surfSlider->setValue(static_cast<int>(_val));
    if(curSurf){
        curSurf->isoval = val;
    }
}

void Data3DWidget::on_surfBut_toggled(bool checked)
{
    if(curSurf){
        curSurf->display = checked;
        if(checked){
            master->addExtraData(&curSurf->gpu_data);
        }else{
            master->delExtraData(&curSurf->gpu_data);
        }
    }else if(checked){
        auto tmp = surfaces.emplace(curData, IsoSurf{
            true,
            static_cast<bool>(ui->surfToggle->checkState()),
            ui->surfVal->text().toDouble(),
            GUI::MeshData{master->getGLGlobals(),
                          {},
                          curData->origin,
                          curData->cell,
                          settings.milCol.val }
            });
        curSurf = &tmp.first->second;
        master->addExtraData(&curSurf->gpu_data);
    }
}
