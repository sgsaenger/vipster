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
    ui->sliceVal->setMaximum(static_cast<int>(curData->extent[
                                              static_cast<size_t>(dir)]));
    ui->surfToggle->setCheckState(Qt::CheckState::Unchecked);
    auto minmax = std::minmax_element(curData->begin(), curData->end());
    auto min = static_cast<double>(*minmax.first);
    auto max = static_cast<double>(*minmax.second);
    auto precision = std::abs(max - min)/1000.;
    auto dec = -static_cast<int>(std::round(std::log10(precision)))+2;
    validator.setRange(min, max, dec);
    ui->maxLabel->setText(QString::number(max));
    ui->minLabel->setText(QString::number(min));
}

void Data3DWidget::on_sliceDir_currentIndexChanged(int index)
{
    ui->sliceVal->setValue(0);
    ui->sliceVal->setMaximum(static_cast<int>(curData->extent[
                                              static_cast<size_t>(index)]));
}

void Data3DWidget::on_sliceVal_valueChanged(int)
{
    // TODO: recreate slice
}

void Data3DWidget::on_surfToggle_stateChanged(int state)
{
    display_pm = state;
    // TODO: reeval isosurface when displayed
}

void Data3DWidget::on_horizontalSlider_valueChanged(int)
{
    //TODO: recreate isosurface, fill surfVal
}

void Data3DWidget::on_surfVal_editingFinished()
{
    //TODO: recreate isosurface, reposition slider
}

//TODO: display something...
void Data3DWidget::on_sliceBut_clicked()
{
}

void Data3DWidget::on_surfBut_clicked()
{
}
