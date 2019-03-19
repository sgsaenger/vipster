#include "data3dwidget.h"
#include "ui_data3dwidget.h"
#include "mainwindow.h"

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

void Data3DWidget::updateWidget(guiChange_t change)
{
    if(change & GuiChange::settings){
        for(auto& p: surfaces){
            p.second.gpu_data.update({{settings.posCol.val,
                                       settings.negCol.val}, 2, 1});
        }
    }
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
        auto isoval = static_cast<double>(curSurf->isoval);
        ui->surfVal->setText(QString::number(isoval));
        auto _val = (isoval - validator.bottom()) *
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

std::vector<GUI::MeshData::Face> mkSlice(size_t dir, float off)
{
    switch(dir){
    case 0:
        return {{{off,0,0},{},{0,0}},{{off,1,0},{},{1,0}},{{off,1,1},{},{1,1}},
                {{off,0,0},{},{0,0}},{{off,0,1},{},{0,1}},{{off,1,1},{},{1,1}}};
    case 1:
        return {{{0,off,0},{},{0,0}},{{1,off,0},{},{1,0}},{{1,off,1},{},{1,1}},
                {{0,off,0},{},{0,0}},{{0,off,1},{},{0,1}},{{1,off,1},{},{1,1}}};
    case 2:
        return {{{0,0,off},{},{0,0}},{{1,0,off},{},{1,0}},{{1,1,off},{},{1,1}},
                {{0,0,off},{},{0,0}},{{0,1,off},{},{0,1}},{{1,1,off},{},{1,1}}};
    default:
        throw Error("Data3DWidget: Invalid direction for data slicing");
    }
}

GUI::MeshData::Texture mkSliceTex(const DataGrid3D_f& dat, size_t dir, size_t off)
{
    GUI::MeshData::Texture texture;
    size_t xl{0}, xh{dat.extent[0]};
    size_t yl{0}, yh{dat.extent[1]};
    size_t zl{0}, zh{dat.extent[2]};
    switch(dir){
    case 0:
        xl = off;
        xh = xl+1;
        texture.width = static_cast<int>(yh);
        texture.height = static_cast<int>(zh);
        break;
    case 1:
        yl = off;
        yh = yl+1;
        texture.width = static_cast<int>(xh);
        texture.height = static_cast<int>(zh);
        break;
    case 2:
        zl = off;
        zh = zl+1;
        texture.width = static_cast<int>(xh);
        texture.height = static_cast<int>(yh);
        break;
    default:
        throw Error("Data3DWidget: Invalid direction for data slicing");
    }
    auto minmax = std::minmax_element(dat.begin(), dat.end());
    auto min = *minmax.first;
    auto max = *minmax.second;
    auto factor = 100/(max-min);
    for(auto z=zl; z < zh; ++z){
        for(auto y=yl; y < yh; ++y){
            for(auto x=xl; x < xh; ++x){
                const auto& val = dat(x,y,z);
                auto tmp = (val-min)*factor;
                texture.data.push_back({static_cast<uint8_t>(std::round(2.55f*tmp)),
                                        static_cast<uint8_t>(std::round(2.55f*(100-abs(2*tmp-100)))),
                                        static_cast<uint8_t>(std::round(2.55f*(100-tmp))),
                                        128});
            }
        }
    }
    return texture;
}

void Data3DWidget::on_sliceDir_currentIndexChanged(int index)
{
    auto _index = static_cast<size_t>(index);
    //block sliceVal from triggering when sliceVal is higher than new max
    QSignalBlocker block{ui->sliceVal};
    ui->sliceVal->setMaximum(static_cast<int>(curData->extent[_index]));
    if(curSlice){
        curSlice->dir = _index;
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
        auto _dir = static_cast<size_t>(ui->sliceDir->currentIndex());
        curSlice->gpu_data.update(mkSlice(_dir, off));
        curSlice->gpu_data.update(mkSliceTex(*curData, _dir, static_cast<size_t>(pos)));
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
                          mkSlice(dir, off),
                          curData->origin,
                          curData->cell,
                          mkSliceTex(*curData, dir, pos)}
            });
        curSlice = &tmp.first->second;
        master->addExtraData(&curSlice->gpu_data);
    }
}

static constexpr Vec vert_off[8]={
    {{0.f, 0.f, 0.f}},
    {{0.f, 0.f, 1.f}},
    {{0.f, 1.f, 0.f}},
    {{0.f, 1.f, 1.f}},
    {{1.f, 0.f, 0.f}},
    {{1.f, 0.f, 1.f}},
    {{1.f, 1.f, 0.f}},
    {{1.f, 1.f, 1.f}},
};

static constexpr int nvert_lut[256]={
    0, 1, 1, 2, 1, 2, 2, 3,
    1, 2, 2, 3, 2, 3, 3, 2,
    1, 2, 2, 3, 2, 3, 3, 4,
    2, 3, 3, 4, 3, 4, 4, 3,
    1, 2, 2, 3, 2, 3, 3, 4,
    2, 3, 3, 4, 3, 4, 4, 3,
    2, 3, 3, 2, 3, 4, 4, 3,
    3, 4, 4, 3, 4, 3, 3, 2,
    1, 2, 2, 3, 2, 3, 3, 4,
    2, 3, 3, 4, 3, 4, 4, 3,
    2, 3, 3, 4, 3, 2, 4, 3,
    3, 4, 4, 3, 4, 3, 3, 2,
    2, 3, 3, 4, 3, 4, 4, 3,
    3, 4, 4, 3, 4, 3, 3, 2,
    3, 4, 4, 3, 4, 3, 3, 2,
    4, 3, 3, 2, 3, 2, 2, 1,
    1, 2, 2, 3, 2, 3, 3, 4,
    2, 3, 3, 4, 3, 4, 4, 3,
    2, 3, 3, 4, 3, 4, 4, 3,
    3, 4, 4, 3, 4, 3, 3, 2,
    2, 3, 3, 4, 3, 4, 4, 3,
    3, 4, 2, 3, 4, 3, 3, 2,
    3, 4, 4, 3, 4, 3, 3, 2,
    4, 3, 3, 2, 3, 2, 2, 1,
    2, 3, 3, 4, 3, 4, 4, 3,
    3, 4, 4, 3, 2, 3, 3, 2,
    3, 4, 4, 3, 4, 3, 3, 2,
    4, 3, 3, 2, 3, 2, 2, 1,
    3, 4, 4, 3, 4, 3, 3, 2,
    4, 3, 3, 2, 3, 2, 2, 1,
    2, 3, 3, 2, 3, 2, 2, 1,
    3, 2, 2, 1, 2, 1, 1, 0};

static constexpr int edge_lut[256]={
    0x000,  0x089,  0x013,  0x09A,  0x980,  0x909,  0x993,  0x91A,
    0x310,  0x399,  0x303,  0x38A,  0xA90,  0xA19,  0xA83,  0xA0A,
    0x04C,  0x0C5,  0x05F,  0x0D6,  0x9CC,  0x945,  0x9DF,  0x956,
    0x35C,  0x3D5,  0x34F,  0x3C6,  0xADC,  0xA55,  0xACF,  0xA46,
    0x026,  0x0AF,  0x035,  0x0BC,  0x9A6,  0x92F,  0x9B5,  0x93C,
    0x336,  0x3BF,  0x325,  0x3AC,  0xAB6,  0xA3F,  0xAA5,  0xA2C,
    0x06A,  0x0E3,  0x079,  0x0F0,  0x9EA,  0x963,  0x9F9,  0x970,
    0x37A,  0x3F3,  0x369,  0x3E0,  0xAFA,  0xA73,  0xAE9,  0xA60,
    0xC40,  0xCC9,  0xC53,  0xCDA,  0x5C0,  0x549,  0x5D3,  0x55A,
    0xF50,  0xFD9,  0xF43,  0xFCA,  0x6D0,  0x659,  0x6C3,  0x64A,
    0xC0C,  0xC85,  0xC1F,  0xC96,  0x58C,  0x505,  0x59F,  0x516,
    0xF1C,  0xF95,  0xF0F,  0xF86,  0x69C,  0x615,  0x68F,  0x606,
    0xC66,  0xCEF,  0xC75,  0xCFC,  0x5E6,  0x56F,  0x5F5,  0x57C,
    0xF76,  0xFFF,  0xF65,  0xFEC,  0x6F6,  0x67F,  0x6E5,  0x66C,
    0xC2A,  0xCA3,  0xC39,  0xCB0,  0x5AA,  0x523,  0x5B9,  0x530,
    0xF3A,  0xFB3,  0xF29,  0xFA0,  0x6BA,  0x633,  0x6A9,  0x620,
    0x620,  0x6A9,  0x633,  0x6BA,  0xFA0,  0xF29,  0xFB3,  0xF3A,
    0x530,  0x5B9,  0x523,  0x5AA,  0xCB0,  0xC39,  0xCA3,  0xC2A,
    0x66C,  0x6E5,  0x67F,  0x6F6,  0xFEC,  0xF65,  0xFFF,  0xF76,
    0x57C,  0x5F5,  0x56F,  0x5E6,  0xCFC,  0xC75,  0xCEF,  0xC66,
    0x606,  0x68F,  0x615,  0x69C,  0xF86,  0xF0F,  0xF95,  0xF1C,
    0x516,  0x59F,  0x505,  0x58C,  0xC96,  0xC1F,  0xC85,  0xC0C,
    0x64A,  0x6C3,  0x659,  0x6D0,  0xFCA,  0xF43,  0xFD9,  0xF50,
    0x65A,  0x5D3,  0x549,  0x5C0,  0xCDA,  0xC53,  0xCC9,  0xC40,
    0xA60,  0xAE9,  0xA73,  0xAFA,  0x3E0,  0x369,  0x3F3,  0x37A,
    0x970,  0x9F9,  0x963,  0x9EA,  0x0F0,  0x079,  0x0E3,  0x06A,
    0xA2C,  0xAA5,  0xA3F,  0xAB6,  0x3AC,  0x325,  0x3BF,  0x336,
    0x93C,  0x9B5,  0x92F,  0x9A6,  0x0BC,  0x035,  0x0AF,  0x026,
    0xA46,  0xACF,  0xA55,  0xADC,  0x3C6,  0x34F,  0x3D5,  0x35C,
    0x956,  0x9DF,  0x945,  0x9CC,  0x0D6,  0x05F,  0x0C5,  0x04C,
    0xA0A,  0xA83,  0xA19,  0xA90,  0x38A,  0x303,  0x399,  0x310,
    0x91A,  0x993,  0x909,  0x980,  0x09A,  0x013,  0x089,  0x000};

static constexpr int tri_lut[256][4][3]={
    {{0, 0, 0},{ 0, 0, 0},{ 0, 0, 0},{ 0, 0, 0}},
    {{0, 3, 7},{ 0, 0, 0},{ 0, 0, 0},{ 0, 0, 0}},
    {{0, 1, 4},{ 0, 0, 0},{ 0, 0, 0},{ 0, 0, 0}},
    {{3, 7, 4},{ 3, 4, 1},{ 0, 0, 0},{ 0, 0, 0}},
    {{7,11, 8},{ 0, 0, 0},{ 0, 0, 0},{ 0, 0, 0}},
    {{0, 3,11},{ 0,11, 8},{ 0, 0, 0},{ 0, 0, 0}},
    {{7,11, 8},{ 0, 1, 4},{ 0, 0, 0},{ 0, 0, 0}},
    {{4, 8, 1},{ 8,11, 1},{ 1,11, 3},{ 0, 0, 0}},
    {{4, 8, 9},{ 0, 0, 0},{ 0, 0, 0},{ 0, 0, 0}},
    {{4, 8, 9},{ 0, 3, 7},{ 0, 0, 0},{ 0, 0, 0}},
    {{0, 1, 8},{ 1, 8, 9},{ 0, 0, 0},{ 0, 0, 0}},
    {{1, 3, 9},{ 3, 9, 8},{ 3, 7, 8},{ 0, 0, 0}},
    {{4, 7,11},{ 4, 9,11},{ 0, 0, 0},{ 0, 0, 0}},
    {{9, 3,11},{ 0, 3, 9},{ 0, 9, 4},{ 0, 0, 0}},
    {{1, 9,11},{ 1, 0,11},{ 0, 7,11},{ 0, 0, 0}},
    {{1, 9,11},{ 1, 3,11},{ 0, 0, 0},{ 0, 0, 0}},
    {{2, 3, 6},{ 0, 0, 0},{ 0, 0, 0},{ 0, 0, 0}},
    {{0, 2, 7},{ 6, 2, 7},{ 0, 0, 0},{ 0, 0, 0}},
    {{2, 3, 6},{ 0, 1, 4},{ 0, 0, 0},{ 0, 0, 0}},
    {{4, 7, 6},{ 4, 6, 1},{ 1, 6, 2},{ 0, 0, 0}},
    {{2, 3, 6},{ 7, 8,11},{ 0, 0, 0},{ 0, 0, 0}},
    {{0, 2, 8},{ 2, 8, 6},{ 8, 6,11},{ 0, 0, 0}},
    {{2, 3, 6},{ 7,11, 8},{ 0, 1, 4},{ 0, 0, 0}},
    {{1, 2, 4},{ 2, 4, 6},{ 4, 6, 8},{ 6, 8,11}},
    {{2, 3, 6},{ 4, 8, 9},{ 0, 0, 0},{ 0, 0, 0}},
    {{0, 2, 7},{ 6, 2, 7},{ 4, 8, 9},{ 0, 0, 0}},
    {{0, 1, 8},{ 1, 8, 9},{ 2, 3, 6},{ 0, 0, 0}},
    {{2, 6, 7},{ 2, 7, 9},{ 1, 2, 9},{ 7, 8, 9}},
    {{2, 3, 6},{ 4, 7,11},{ 4, 9,11},{ 0, 0, 0}},
    {{2, 6, 0},{ 0, 6, 9},{ 0, 9, 4},{ 6, 9,11}},
    {{1, 9,11},{ 1, 0,11},{ 0, 7,11},{ 2, 3, 6}},
    {{1, 9,11},{ 1, 2,11},{ 2, 6,11},{ 0, 0, 0}},
    {{1, 2, 5},{ 0, 0, 0},{ 0, 0, 0},{ 0, 0, 0}},
    {{0, 3, 7},{ 1, 2, 5},{ 0, 0, 0},{ 0, 0, 0}},
    {{0, 2, 5},{ 0, 4, 5},{ 0, 0, 0},{ 0, 0, 0}},
    {{4, 5, 7},{ 5, 7, 2},{ 2, 7, 3},{ 0, 0, 0}},
    {{1, 2, 5},{ 7,11, 8},{ 0, 0, 0},{ 0, 0, 0}},
    {{0, 3,11},{ 0,11, 8},{ 1, 2, 5},{ 0, 0, 0}},
    {{0, 2, 5},{ 0, 4, 5},{ 7,11, 8},{ 0, 0, 0}},
    {{2, 3,11},{ 2,11, 4},{ 2, 4, 5},{ 4, 8,11}},
    {{1, 2, 5},{ 4, 8, 9},{ 0, 0, 0},{ 0, 0, 0}},
    {{0, 3, 7},{ 1, 2, 5},{ 4, 8, 9},{ 0, 0, 0}},
    {{0, 2, 8},{ 2, 8, 5},{ 5, 8, 9},{ 0, 0, 0}},
    {{3, 7, 2},{ 2, 7, 5},{ 5, 7, 8},{ 5, 8, 9}},
    {{4, 7,11},{ 4, 9,11},{ 1, 2, 5},{ 0, 0, 0}},
    {{1, 2, 5},{ 3, 9,11},{ 3, 9, 4},{ 3, 4, 0}},
    {{0, 2, 5},{ 0, 5,11},{ 5,11, 9},{ 0, 7,11}},
    {{3, 9,11},{ 3, 9, 5},{ 3, 5, 2},{ 0, 0, 0}},
    {{1, 3, 6},{ 1, 5, 6},{ 0, 0, 0},{ 0, 0, 0}},
    {{5, 6, 7},{ 5, 7, 0},{ 5, 0, 1},{ 0, 0, 0}},
    {{4, 5, 6},{ 4, 6, 3},{ 4, 3, 0},{ 0, 0, 0}},
    {{4, 5, 6},{ 4, 6, 7},{ 0, 0, 0},{ 0, 0, 0}},
    {{7, 8,11},{ 1, 3, 5},{ 3, 5, 6},{ 0, 0, 0}},
    {{1, 5, 6},{ 1, 6, 8},{ 1, 8, 0},{ 6, 8,11}},
    {{4, 5, 6},{ 4, 6, 3},{ 4, 3, 0},{ 7, 8,11}},
    {{4, 5, 6},{ 4, 6,11},{ 4, 8,11},{ 0, 0, 0}},
    {{4, 8, 9},{ 1, 3, 6},{ 1, 5, 6},{ 0, 0, 0}},
    {{5, 6, 7},{ 5, 7, 0},{ 5, 0, 1},{ 4, 8, 9}},
    {{0, 8, 9},{ 0, 6, 9},{ 9, 5, 6},{ 0, 6, 3}},
    {{5, 6, 7},{ 5, 7, 8},{ 5, 8, 9},{ 0, 0, 0}},
    {{1, 3, 5},{ 3, 5, 6},{ 4, 7, 9},{ 7, 9,11}},
    {{0, 1, 4},{ 5, 6, 9},{ 6, 9,11},{ 0, 0, 0}},
    {{0, 3, 7},{ 5, 6, 9},{ 6, 9,11},{ 0, 0, 0}},
    {{5, 6, 9},{ 6, 9,11},{ 0, 0, 0},{ 0, 0, 0}},
    {{6,10,11},{ 0, 0, 0},{ 0, 0, 0},{ 0, 0, 0}},
    {{0, 3, 7},{ 6,10,11},{ 0, 0, 0},{ 0, 0, 0}},
    {{0, 1, 4},{ 6,10,11},{ 0, 0, 0},{ 0, 0, 0}},
    {{1, 3, 4},{ 3, 4, 7},{ 6,10,11},{ 0, 0, 0}},
    {{6, 7, 8},{ 6, 8,10},{ 0, 0, 0},{ 0, 0, 0}},
    {{0, 8,10},{ 0,10, 6},{ 0, 6, 3},{ 0, 0, 0}},
    {{6, 7, 8},{ 6, 8,10},{ 0, 1, 4},{ 0, 0, 0}},
    {{1, 3, 4},{ 3, 4,10},{ 3,10, 6},{ 4, 8,10}},
    {{4, 8, 9},{ 6,10,11},{ 0, 0, 0},{ 0, 0, 0}},
    {{0, 3, 7},{ 4, 8, 9},{ 6,10,11},{ 0, 0, 0}},
    {{0, 1, 8},{ 1, 8, 9},{ 6,10,11},{ 0, 0, 0}},
    {{1, 3, 9},{ 3, 9, 8},{ 3, 8, 7},{ 6,10,11}},
    {{4, 6, 7},{ 4, 6,10},{ 4,10, 9},{ 0, 0, 0}},
    {{0, 3, 4},{ 3, 4, 9},{ 3, 9, 6},{ 6, 9,10}},
    {{0, 1, 9},{ 0, 9, 6},{ 0, 6, 7},{ 6, 9,10}},
    {{1, 3, 9},{ 3, 9, 6},{ 6, 9,10},{ 0, 0, 0}},
    {{2, 3,10},{ 3,10,11},{ 0, 0, 0},{ 0, 0, 0}},
    {{0, 2,10},{ 0,10, 7},{ 7,10,11},{ 0, 0, 0}},
    {{2, 3,10},{ 3,10,11},{ 0, 1, 4},{ 0, 0, 0}},
    {{1, 4, 7},{ 1, 7,10},{ 7,10,11},{ 1,10, 2}},
    {{2, 8,10},{ 2, 8, 7},{ 2, 7, 3},{ 0, 0, 0}},
    {{0, 2, 8},{ 2, 8,10},{ 0, 0, 0},{ 0, 0, 0}},
    {{2, 8,10},{ 2, 8, 7},{ 2, 7, 3},{ 0, 1, 4}},
    {{2, 8,10},{ 2, 8, 4},{ 2, 4, 1},{ 0, 0, 0}},
    {{2, 3,10},{ 3,10,11},{ 4, 8, 9},{ 0, 0, 0}},
    {{0, 2,10},{ 0,10,11},{ 0,11, 7},{ 4, 8, 9}},
    {{2, 3,10},{ 3,10,11},{ 0, 1, 8},{ 1, 8, 9}},
    {{1, 2, 9},{ 2, 9,10},{ 7, 8,11},{ 0, 0, 0}},
    {{4, 9, 7},{ 7, 9, 2},{ 2, 9,10},{ 2, 3, 7}},
    {{0, 2,10},{ 0,10, 9},{ 0, 9, 4},{ 0, 0, 0}},
    {{0, 3, 7},{ 1, 2, 9},{ 2, 9,10},{ 0, 0, 0}},
    {{1, 2, 9},{ 2, 9,10},{ 0, 0, 0},{ 0, 0, 0}},
    {{1, 2, 5},{ 6,10,11},{ 0, 0, 0},{ 0, 0, 0}},
    {{1, 2, 5},{ 6,10,11},{ 0, 3, 7},{ 0, 0, 0}},
    {{0, 4, 5},{ 0, 5, 2},{ 6,10,11},{ 0, 0, 0}},
    {{4, 5, 7},{ 5, 7, 3},{ 5, 3, 2},{ 6,10,11}},
    {{1, 2, 5},{ 6, 7, 8},{ 6, 8,10},{ 0, 0, 0}},
    {{1, 2, 5},{ 0, 8,10},{ 0,10, 3},{ 3,10, 6}},
    {{0, 4, 5},{ 0, 5, 2},{ 6, 7, 8},{ 6, 8,10}},
    {{2, 3, 6},{ 4, 5, 8},{ 5, 8,10},{ 0, 0, 0}},
    {{1, 2, 5},{ 4, 8, 9},{ 6,10,11},{ 0, 0, 0}},
    {{0, 3, 7},{ 1, 2, 5},{ 4, 8, 9},{ 6,10,11}},
    {{0, 2, 8},{ 2, 8, 5},{ 5, 8, 9},{ 6,10,11}},
    {{2, 3, 6},{ 7, 8,11},{ 5, 9,10},{ 0, 0, 0}},
    {{4, 6, 7},{ 4, 6, 9},{ 6, 9,10},{ 1, 2, 5}},
    {{0, 1, 4},{ 5, 9,10},{ 2, 3, 6},{ 0, 0, 0}},
    {{0, 2, 7},{ 2, 7, 6},{ 5, 9,10},{ 0, 0, 0}},
    {{2, 3, 6},{ 5, 9,10},{ 0, 0, 0},{ 0, 0, 0}},
    {{1, 3,11},{ 1,11,10},{ 1,10, 5},{ 0, 0, 0}},
    {{0, 1, 7},{ 1, 7, 5},{ 7, 5,11},{ 5,10,11}},
    {{0, 4, 5},{ 0, 5,11},{ 0, 3,11},{ 5,11,10}},
    {{4, 5, 7},{ 5, 7,10},{ 7,10,11},{ 0, 0, 0}},
    {{1, 5, 3},{ 3, 5, 8},{ 3, 8, 7},{ 5, 8,10}},
    {{0, 8,10},{ 0,10, 5},{ 0, 5, 1},{ 0, 0, 0}},
    {{0, 3, 7},{ 4, 5, 8},{ 5, 8,10},{ 0, 0, 0}},
    {{4, 5, 8},{ 5, 8,10},{ 0, 0, 0},{ 0, 0, 0}},
    {{4, 8, 9},{ 1, 3,11},{ 1,11,10},{ 1,10, 5}},
    {{0, 1, 4},{ 7, 8,11},{ 5, 9,10},{ 0, 0, 0}},
    {{0, 3, 8},{ 3, 8,11},{ 5, 9,10},{ 0, 0, 0}},
    {{7, 8,11},{ 5, 9,10},{ 0, 0, 0},{ 0, 0, 0}},
    {{1, 3, 4},{ 3, 4, 7},{ 5, 9,10},{ 0, 0, 0}},
    {{0, 1, 4},{ 5, 9,10},{ 0, 0, 0},{ 0, 0, 0}},
    {{0, 3, 7},{ 5, 9,10},{ 0, 0, 0},{ 0, 0, 0}},
    {{5, 9,10},{ 0, 0, 0},{ 0, 0, 0},{ 0, 0, 0}},
    {{5, 9,10},{ 0, 0, 0},{ 0, 0, 0},{ 0, 0, 0}},
    {{0, 3, 7},{ 5, 9,10},{ 0, 0, 0},{ 0, 0, 0}},
    {{0, 1, 4},{ 5, 9,10},{ 0, 0, 0},{ 0, 0, 0}},
    {{1, 3, 4},{ 3, 4, 7},{ 5, 9,10},{ 0, 0, 0}},
    {{7, 8,11},{ 5, 9,10},{ 0, 0, 0},{ 0, 0, 0}},
    {{0, 3, 8},{ 3, 8,11},{ 5, 9,10},{ 0, 0, 0}},
    {{0, 1, 4},{ 7, 8,11},{ 5, 9,10},{ 0, 0, 0}},
    {{4, 8, 9},{ 1, 3,11},{ 1,11,10},{ 1,10, 5}},
    {{4, 5, 8},{ 5, 8,10},{ 0, 0, 0},{ 0, 0, 0}},
    {{0, 3, 7},{ 4, 5, 8},{ 5, 8,10},{ 0, 0, 0}},
    {{0, 8,10},{ 0,10, 5},{ 0, 5, 1},{ 0, 0, 0}},
    {{1, 5, 3},{ 3, 5, 8},{ 3, 8, 7},{ 5, 8,10}},
    {{4, 5, 7},{ 5, 7,10},{ 7,10,11},{ 0, 0, 0}},
    {{0, 4, 5},{ 0, 5,11},{ 0, 3,11},{ 5,11,10}},
    {{0, 1, 7},{ 1, 7, 5},{ 7, 5,11},{ 5,10,11}},
    {{1, 3,11},{ 1,11,10},{ 1,10, 5},{ 0, 0, 0}},
    {{2, 3, 6},{ 5, 9,10},{ 0, 0, 0},{ 0, 0, 0}},
    {{0, 2, 7},{ 2, 7, 6},{ 5, 9,10},{ 0, 0, 0}},
    {{0, 1, 4},{ 5, 9,10},{ 2, 3, 6},{ 0, 0, 0}},
    {{4, 6, 7},{ 4, 6, 9},{ 6, 9,10},{ 1, 2, 5}},
    {{2, 3, 6},{ 7, 8,11},{ 5, 9,10},{ 0, 0, 0}},
    {{0, 2, 8},{ 2, 8, 5},{ 5, 8, 9},{ 6,10,11}},
    {{0, 3, 7},{ 1, 2, 5},{ 4, 8, 9},{ 6,10,11}},
    {{1, 2, 5},{ 4, 8, 9},{ 6,10,11},{ 0, 0, 0}},
    {{2, 3, 6},{ 4, 5, 8},{ 5, 8,10},{ 0, 0, 0}},
    {{0, 4, 5},{ 0, 5, 2},{ 6, 7, 8},{ 6, 8,10}},
    {{1, 2, 5},{ 0, 8,10},{ 0,10, 3},{ 3,10, 6}},
    {{1, 2, 5},{ 6, 7, 8},{ 6, 8,10},{ 0, 0, 0}},
    {{4, 5, 7},{ 5, 7, 3},{ 5, 3, 2},{ 6,10,11}},
    {{0, 4, 5},{ 0, 5, 2},{ 6,10,11},{ 0, 0, 0}},
    {{1, 2, 5},{ 6,10,11},{ 0, 3, 7},{ 0, 0, 0}},
    {{1, 2, 5},{ 6,10,11},{ 0, 0, 0},{ 0, 0, 0}},
    {{1, 2, 9},{ 2, 9,10},{ 0, 0, 0},{ 0, 0, 0}},
    {{0, 3, 7},{ 1, 2, 9},{ 2, 9,10},{ 0, 0, 0}},
    {{0, 2,10},{ 0,10, 9},{ 0, 9, 4},{ 0, 0, 0}},
    {{4, 9, 7},{ 7, 9, 2},{ 2, 9,10},{ 2, 3, 7}},
    {{1, 2, 9},{ 2, 9,10},{ 7, 8,11},{ 0, 0, 0}},
    {{2, 3,10},{ 3,10,11},{ 0, 1, 8},{ 1, 8, 9}},
    {{0, 2,10},{ 0,10,11},{ 0,11, 7},{ 4, 8, 9}},
    {{2, 3,10},{ 3,10,11},{ 4, 8, 9},{ 0, 0, 0}},
    {{2, 8,10},{ 2, 8, 4},{ 2, 4, 1},{ 0, 0, 0}},
    {{2, 8,10},{ 2, 8, 7},{ 2, 7, 3},{ 0, 1, 4}},
    {{0, 2, 8},{ 2, 8,10},{ 0, 0, 0},{ 0, 0, 0}},
    {{2, 8,10},{ 2, 8, 7},{ 2, 7, 3},{ 0, 0, 0}},
    {{1, 4, 7},{ 1, 7,10},{ 7,10,11},{ 1,10, 2}},
    {{2, 3,10},{ 3,10,11},{ 0, 1, 4},{ 0, 0, 0}},
    {{0, 2,10},{ 0,10, 7},{ 7,10,11},{ 0, 0, 0}},
    {{2, 3,10},{ 3,10,11},{ 0, 0, 0},{ 0, 0, 0}},
    {{1, 3, 9},{ 3, 9, 6},{ 6, 9,10},{ 0, 0, 0}},
    {{0, 1, 9},{ 0, 9, 6},{ 0, 6, 7},{ 6, 9,10}},
    {{0, 3, 4},{ 3, 4, 9},{ 3, 9, 6},{ 6, 9,10}},
    {{4, 6, 7},{ 4, 6,10},{ 4,10, 9},{ 0, 0, 0}},
    {{1, 3, 9},{ 3, 9, 8},{ 3, 8, 7},{ 6,10,11}},
    {{0, 1, 8},{ 1, 8, 9},{ 6,10,11},{ 0, 0, 0}},
    {{0, 3, 7},{ 4, 8, 9},{ 6,10,11},{ 0, 0, 0}},
    {{4, 8, 9},{ 6,10,11},{ 0, 0, 0},{ 0, 0, 0}},
    {{1, 3, 4},{ 3, 4,10},{ 3,10, 6},{ 4, 9,10}},
    {{6, 7, 8},{ 6, 8,10},{ 0, 1, 4},{ 0, 0, 0}},
    {{0, 8,10},{ 0,10, 6},{ 0, 6, 3},{ 0, 0, 0}},
    {{6, 7, 8},{ 6, 8,10},{ 0, 0, 0},{ 0, 0, 0}},
    {{1, 3, 4},{ 3, 4, 7},{ 6,10,11},{ 0, 0, 0}},
    {{0, 1, 4},{ 6,10,11},{ 0, 0, 0},{ 0, 0, 0}},
    {{0, 3, 7},{ 6,10,11},{ 0, 0, 0},{ 0, 0, 0}},
    {{6,10,11},{ 0, 0, 0},{ 0, 0, 0},{ 0, 0, 0}},
    {{5, 6, 9},{ 6, 9,11},{ 0, 0, 0},{ 0, 0, 0}},
    {{0, 3, 7},{ 5, 6, 9},{ 6, 9,11},{ 0, 0, 0}},
    {{0, 1, 4},{ 5, 6, 9},{ 6, 9,11},{ 0, 0, 0}},
    {{1, 3, 5},{ 3, 5, 6},{ 4, 7, 9},{ 7, 9,11}},
    {{5, 6, 7},{ 5, 7, 8},{ 5, 8, 9},{ 0, 0, 0}},
    {{0, 8, 9},{ 0, 6, 9},{ 9, 5, 6},{ 0, 6, 3}},
    {{5, 6, 7},{ 5, 7, 0},{ 5, 0, 1},{ 4, 8, 9}},
    {{4, 8, 9},{ 1, 3, 6},{ 1, 5, 6},{ 0, 0, 0}},
    {{4, 5, 6},{ 4, 6,11},{ 4, 8,11},{ 0, 0, 0}},
    {{4, 5, 6},{ 4, 6, 3},{ 4, 3, 0},{ 7, 8,11}},
    {{1, 5, 6},{ 1, 6, 8},{ 1, 8, 0},{ 6, 8,11}},
    {{7, 8,11},{ 1, 3, 5},{ 3, 5, 6},{ 0, 0, 0}},
    {{4, 5, 6},{ 4, 6, 7},{ 0, 0, 0},{ 0, 0, 0}},
    {{4, 5, 6},{ 4, 6, 3},{ 4, 3, 0},{ 0, 0, 0}},
    {{5, 6, 7},{ 5, 7, 0},{ 5, 0, 1},{ 0, 0, 0}},
    {{1, 3, 6},{ 1, 5, 6},{ 0, 0, 0},{ 0, 0, 0}},
    {{3, 9,11},{ 3, 9, 5},{ 3, 5, 2},{ 0, 0, 0}},
    {{0, 2, 5},{ 0, 5,11},{ 5,11, 9},{ 0, 7,11}},
    {{1, 2, 5},{ 3, 9,11},{ 3, 9, 4},{ 3, 4, 0}},
    {{4, 7,11},{ 4, 9,11},{ 1, 2, 5},{ 0, 0, 0}},
    {{3, 7, 2},{ 2, 7, 5},{ 5, 7, 8},{ 5, 8, 9}},
    {{0, 2, 8},{ 2, 8, 5},{ 5, 8, 9},{ 0, 0, 0}},
    {{0, 3, 7},{ 1, 2, 5},{ 4, 8, 9},{ 0, 0, 0}},
    {{1, 2, 5},{ 4, 8, 9},{ 0, 0, 0},{ 0, 0, 0}},
    {{2, 3,11},{ 2,11, 4},{ 2, 4, 5},{ 4, 8,11}},
    {{0, 2, 5},{ 0, 4, 5},{ 7,11, 8},{ 0, 0, 0}},
    {{0, 3,11},{ 0,11, 8},{ 1, 2, 5},{ 0, 0, 0}},
    {{1, 2, 5},{ 7,11, 8},{ 0, 0, 0},{ 0, 0, 0}},
    {{4, 5, 7},{ 5, 7, 2},{ 2, 7, 3},{ 0, 0, 0}},
    {{0, 2, 5},{ 0, 4, 5},{ 0, 0, 0},{ 0, 0, 0}},
    {{0, 3, 7},{ 1, 2, 5},{ 0, 0, 0},{ 0, 0, 0}},
    {{1, 2, 5},{ 0, 0, 0},{ 0, 0, 0},{ 0, 0, 0}},
    {{1, 9,11},{ 1, 2,11},{ 2, 6,11},{ 0, 0, 0}},
    {{1, 9,11},{ 1, 0,11},{ 0, 7,11},{ 2, 3, 6}},
    {{2, 6, 0},{ 0, 6, 9},{ 0, 9, 4},{ 6, 9,11}},
    {{2, 3, 6},{ 4, 7,11},{ 4, 9,11},{ 0, 0, 0}},
    {{2, 6, 7},{ 2, 7, 9},{ 1, 2, 9},{ 7, 8, 9}},
    {{0, 1, 8},{ 1, 8, 9},{ 2, 3, 6},{ 0, 0, 0}},
    {{0, 2, 7},{ 6, 2, 7},{ 4, 8, 9},{ 0, 0, 0}},
    {{2, 3, 6},{ 4, 8, 9},{ 0, 0, 0},{ 0, 0, 0}},
    {{1, 2, 4},{ 2, 4, 6},{ 4, 6, 8},{ 6, 8,11}},
    {{2, 3, 6},{ 7,11, 8},{ 0, 1, 4},{ 0, 0, 0}},
    {{0, 2, 8},{ 2, 8, 6},{ 8, 6,11},{ 0, 0, 0}},
    {{2, 3, 6},{ 7, 8,11},{ 0, 0, 0},{ 0, 0, 0}},
    {{4, 7, 6},{ 4, 6, 1},{ 1, 6, 2},{ 0, 0, 0}},
    {{2, 3, 6},{ 0, 1, 4},{ 0, 0, 0},{ 0, 0, 0}},
    {{0, 2, 7},{ 6, 2, 7},{ 0, 0, 0},{ 0, 0, 0}},
    {{2, 3, 6},{ 0, 0, 0},{ 0, 0, 0},{ 0, 0, 0}},
    {{1, 9,11},{ 1, 3,11},{ 0, 0, 0},{ 0, 0, 0}},
    {{1, 9,11},{ 1, 0,11},{ 0, 7,11},{ 0, 0, 0}},
    {{9, 3,11},{ 0, 3, 9},{ 0, 9, 4},{ 0, 0, 0}},
    {{4, 7,11},{ 4, 9,11},{ 0, 0, 0},{ 0, 0, 0}},
    {{1, 3, 9},{ 3, 9, 8},{ 3, 7, 8},{ 0, 0, 0}},
    {{0, 1, 8},{ 1, 8, 9},{ 0, 0, 0},{ 0, 0, 0}},
    {{4, 8, 9},{ 0, 3, 7},{ 0, 0, 0},{ 0, 0, 0}},
    {{4, 8, 9},{ 0, 0, 0},{ 0, 0, 0},{ 0, 0, 0}},
    {{4, 8, 1},{ 8,11, 1},{ 1,11, 3},{ 0, 0, 0}},
    {{7,11, 8},{ 0, 1, 4},{ 0, 0, 0},{ 0, 0, 0}},
    {{0, 3,11},{ 0,11, 8},{ 0, 0, 0},{ 0, 0, 0}},
    {{7,11, 8},{ 0, 0, 0},{ 0, 0, 0},{ 0, 0, 0}},
    {{3, 7, 4},{ 3, 4, 1},{ 0, 0, 0},{ 0, 0, 0}},
    {{0, 1, 4},{ 0, 0, 0},{ 0, 0, 0},{ 0, 0, 0}},
    {{0, 3, 7},{ 0, 0, 0},{ 0, 0, 0},{ 0, 0, 0}},
    {{0, 0, 0},{ 0, 0, 0},{ 0, 0, 0},{ 0, 0, 0}}
};

DataGrid3D_v makeGradient(const DataGrid3D_f& dat)
{
    DataGrid3D_v grad{dat.extent};
    auto x = dat.extent[0];
    auto y = dat.extent[1];
    auto z = dat.extent[2];
    size_t il,ih,jl,jh,kl,kh;
    for(size_t i=0; i<x; ++i){
        if(i == 0){
            il = x-1; ih = 1;
        }else if(i == x-1){
            il = x-2; ih = 0;
        }else{
            il = i-1; ih = i+1;
        }
        for(size_t j=0; j<y; ++j){
            if(j == 0){
                jl = y-1; jh = 1;
            }else if(j == y-1){
                jl = y-2; jh = 0;
            }else{
                jl = j-1; jh = j+1;
            }
            for(size_t k=0; k<z; ++k){
                if(k == 0){
                    kl = z-1; kh = 1;
                }else if(k == z-1){
                    kl = z-2; kh = 0;
                }else{
                    kl = k-1; kh = k+1;
                }
                grad(i, j, k) = {
                        dat(il, j , k ) - dat(ih, j , k ),
                        dat(i , jl, k ) - dat(i , jh, k ),
                        dat(i , j , kl) - dat(i , j , kh),
                };
            }
        }
    }
    return grad;
}

std::vector<GUI::MeshData::Face> marchingCubes(const DataGrid3D_f& dat, float isoval)
{
    std::vector<GUI::MeshData::Face> faces;
    auto grad = makeGradient(dat);
    auto x = dat.extent[0];
    auto y = dat.extent[1];
    auto z = dat.extent[2];
    float normdir = isoval<0 ? -1.f: 1.f;
    decltype(GUI::MeshData::Face::uv)::value_type color_u = isoval<0 ? 1 : 0;

    faces.reserve(4*x*y*z);

    // tmp-variables for multiple interpolation-invocations
    Vec tmppos[12], tmpnorm[12];

    auto interpol = [&](size_t n, size_t e1, size_t e2,
                        size_t i_1, size_t j_1, size_t k_1,
                        size_t i_2, size_t j_2, size_t k_2){
        auto ratio = (isoval - dat(i_1, j_1, k_1)) / (dat(i_2, j_2, k_2) - dat(i_1, j_1, k_1));
        tmppos[n] = vert_off[e1] + ratio*(vert_off[e2] - vert_off[e1]);
        tmpnorm[n] = grad(i_1, j_1, k_1) + ratio*(grad(i_2, j_2, k_2) - grad(i_1, j_1, k_1));
    };

    for(size_t i=0; i<x; ++i){
        size_t i2 = i==x-1 ? 0 : i+1;
        for(size_t j=0; j<y; ++j){
            size_t j2 = j==y-1 ? 0 : j+1;
            for(size_t k=0; k<z; ++k){
                size_t k2 = k==z-1 ? 0 : k+1;
                uint8_t vert_sum{0};
                if (dat(i , j , k )<isoval) vert_sum |= 0x01;
                if (dat(i , j , k2)<isoval) vert_sum |= 0x02;
                if (dat(i , j2, k )<isoval) vert_sum |= 0x04;
                if (dat(i , j2, k2)<isoval) vert_sum |= 0x08;
                if (dat(i2, j , k )<isoval) vert_sum |= 0x10;
                if (dat(i2, j , k2)<isoval) vert_sum |= 0x20;
                if (dat(i2, j2, k )<isoval) vert_sum |= 0x40;
                if (dat(i2, j2, k2)<isoval) vert_sum |= 0x80;

                if (edge_lut[vert_sum]&0x001) interpol( 0,0,1, i , j , k , i , j , k2);
                if (edge_lut[vert_sum]&0x002) interpol( 1,1,5, i , j , k2, i2, j , k2);
                if (edge_lut[vert_sum]&0x004) interpol( 2,5,4, i2, j , k2, i2, j , k );
                if (edge_lut[vert_sum]&0x008) interpol( 3,0,4, i , j , k , i2, j , k );
                if (edge_lut[vert_sum]&0x010) interpol( 4,1,3, i , j , k2, i , j2, k2);
                if (edge_lut[vert_sum]&0x020) interpol( 5,5,7, i2, j , k2, i2, j2, k2);
                if (edge_lut[vert_sum]&0x040) interpol( 6,4,6, i2, j , k , i2, j2, k );
                if (edge_lut[vert_sum]&0x080) interpol( 7,0,2, i , j , k , i , j2, k );
                if (edge_lut[vert_sum]&0x100) interpol( 8,2,3, i , j2, k , i , j2, k2);
                if (edge_lut[vert_sum]&0x200) interpol( 9,3,7, i , j2, k2, i2, j2, k2);
                if (edge_lut[vert_sum]&0x400) interpol(10,6,7, i2, j2, k , i2, j2, k2);
                if (edge_lut[vert_sum]&0x800) interpol(11,2,6, i , j2, k , i2, j2, k );

                for(int l=0; l<nvert_lut[vert_sum]; ++l){
                    for(const auto& vert: tri_lut[vert_sum][l]){
                        faces.push_back({Vec{(i + tmppos[vert][0])/x,
                                             (j + tmppos[vert][1])/y,
                                             (k + tmppos[vert][2])/z},
                                         normdir * tmpnorm[vert],
                                         {color_u,0}});
                    }
                }
            }
        }
    }
    return faces;
}

std::vector<GUI::MeshData::Face> Data3DWidget::mkSurf(float isoval, bool pm)
{
    std::vector<GUI::MeshData::Face> retval;
    if(isoval > static_cast<float>(validator.top()) ||
       isoval < static_cast<float>(validator.bottom())){
        retval = {};
    }else{
        retval = marchingCubes(*curData, isoval);
    }
    if(pm){
        auto tmp = mkSurf(-isoval, false);
        retval.insert(retval.end(), tmp.begin(), tmp.end());
    }
    return retval;
}

void Data3DWidget::on_surfToggle_stateChanged(int state)
{
    if(curSurf){
        curSurf->plusmin = state;
        curSurf->gpu_data.update(mkSurf(curSurf->isoval, curSurf->plusmin));
        if(curSurf->display){
            triggerUpdate(GuiChange::extra);
        }
    }
}

void Data3DWidget::on_surfSlider_valueChanged(int val)
{
    QSignalBlocker block{ui->surfVal};
    auto _val = (val * (validator.top()-validator.bottom()) /
                 ui->surfSlider->maximum()) + validator.bottom();
    ui->surfVal->setText(QString::number(_val));
    if(curSurf){
        curSurf->isoval = static_cast<float>(_val);
        curSurf->gpu_data.update(mkSurf(curSurf->isoval, curSurf->plusmin));
        if(curSurf->display){
            triggerUpdate(GuiChange::extra);
        }
    }
}

void Data3DWidget::on_surfVal_editingFinished()
{
    QSignalBlocker block{ui->surfSlider};
    auto val = ui->surfVal->text().toFloat();
    auto _val = (static_cast<double>(val) - validator.bottom()) *
                ui->surfSlider->maximum() /
                (validator.top()-validator.bottom());
    ui->surfSlider->setValue(static_cast<int>(_val));
    if(curSurf){
        curSurf->isoval = val;
        curSurf->gpu_data.update(mkSurf(val, curSurf->plusmin));
        if(curSurf->display){
            triggerUpdate(GuiChange::extra);
        }
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
        auto isoval = ui->surfVal->text().toFloat();
        auto pm = static_cast<bool>(ui->surfToggle->checkState());
        auto tmp = surfaces.emplace(curData, IsoSurf{
            true, pm, isoval,
            GUI::MeshData{master->getGLGlobals(),
                          mkSurf(isoval, pm),
                          curData->origin,
                          curData->cell,
                          {{settings.posCol.val,
                            settings.negCol.val}, 2, 1}}
            });
        curSurf = &tmp.first->second;
        master->addExtraData(&curSurf->gpu_data);
    }
}
