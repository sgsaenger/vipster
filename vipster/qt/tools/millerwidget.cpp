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
    if((change & guiStepChanged) == guiStepChanged){
        curStep = master->curStep;
        auto pos = planes.find(curStep);
        QSignalBlocker block{this};
        if(pos != planes.end()){
            curPlane = &pos->second;
            ui->hSel->setValue(curPlane->hkl[0]);
            ui->kSel->setValue(curPlane->hkl[1]);
            ui->lSel->setValue(curPlane->hkl[2]);
            ui->xOff->setValue(static_cast<double>(curPlane->offset[0]));
            ui->yOff->setValue(static_cast<double>(curPlane->offset[1]));
            ui->zOff->setValue(static_cast<double>(curPlane->offset[2]));
            ui->pushButton->setChecked(curPlane->display);
            active = true;
            if(active){
                //TODO: push data through to guiwrapper!
            }
        }else{
            curPlane = nullptr;
            ui->hSel->setValue(0);
            ui->kSel->setValue(0);
            ui->lSel->setValue(0);
            ui->xOff->setValue(0.);
            ui->yOff->setValue(0.);
            ui->zOff->setValue(0.);
            ui->pushButton->setChecked(false);
            active = false;
        }
    }else if(change & GuiChange::cell){
        // TODO: recalc vertices when needed
    }
}

void MillerWidget::updateIndex(int idx)
{
//    const auto* sender = QObject::sender();
//    if(sender == ui->hSel){
//        plane[0] = static_cast<uint8_t>(idx);
//    }else if(sender == ui->kSel){
//        plane[1] = static_cast<uint8_t>(idx);
//    }else if(sender == ui->lSel){
//        plane[2] = static_cast<uint8_t>(idx);
//    }else{
//        throw Error("Unknown sender for HKL-plane index");
//    }
//    if(display && ((plane[0] == 0) || (plane[1] == 0) || (plane[2] == 0))){
//        display = false;
//        triggerUpdate(GuiChange::extra);
//    }
}

void MillerWidget::updateOffset(double off)
{
//    const auto* sender = QObject::sender();
//    if(sender == ui->xOff){
//        offset[0] = static_cast<float>(off);
//    }else if(sender == ui->yOff){
//        offset[1] = static_cast<float>(off);
//    }else if(sender == ui->zOff){
//        offset[2] = static_cast<float>(off);
//    }else{
//        throw Error("Unknown sender for HKL-plane offset");
//    }
}

void MillerWidget::toggleView()
{
}
