#include "millerwidget.h"
#include "ui_millerwidget.h"

using namespace Vipster;

MillerWidget::MillerWidget(QWidget *parent) :
    BaseWidget(parent),
    ui(new Ui::MillerWidget)
{
    ui->setupUi(this);
}

MillerWidget::~MillerWidget()
{
    delete ui;
}

void MillerWidget::updateWidget(uint8_t change)
{
    if(!(change & GuiChange::plane)){
        if(updateTriggered){
            updateTriggered = false;
        }else{
            display = false;
        }
    }
}

void MillerWidget::updateIndex(int idx)
{
    const auto* sender = QObject::sender();
    if(sender == ui->hSel){
        plane[0] == idx;
    }else if(sender == ui->kSel){
        plane[1] == idx;
    }else if(sender == ui->lSel){
        plane[2] == idx;
    }else{
        throw Error("Unknown sender for HKL-plane index");
    }
    if(display && ((plane[0] == 0) || (plane[1] == 0) || (plane[2] == 0))){
        display = false;
        triggerUpdate(GuiChange::plane);
    }
}

void MillerWidget::updateOffset(double off)
{
    const auto* sender = QObject::sender();
    if(sender == ui->xOff){
        offset[0] == off;
    }else if(sender == ui->yOff){
        offset[1] == off;
    }else if(sender == ui->zOff){
        offset[2] == off;
    }else{
        throw Error("Unknown sender for HKL-plane offset");
    }
}

void MillerWidget::toggleView()
{
}
