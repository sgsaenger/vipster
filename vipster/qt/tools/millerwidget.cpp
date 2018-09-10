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

std::vector<Vec> mkVertices(const StepProper* step,
                            const std::array<int8_t,3>& hkl)
{
    std::vector<Vec> vert{};
    auto hasH = hkl[0] != 0;
    auto hasK = hkl[1] != 0;
    auto hasL = hkl[2] != 0;
    float valX = 1.f/hkl[0];
    float valY = 1.f/hkl[1];
    float valZ = 1.f/hkl[2];
    if(hasH){
        if(hasK){
            if(hasL){
                // hkl
                for(int h=0; h < abs(hkl[0]); ++h){
                    for(int k=0; k < abs(hkl[1]); ++k){
                        for(int l=0; l < abs(hkl[2]); ++l){
                            vert.push_back({(h+1)*valX, k*valY,     l*valZ});
                            vert.push_back({h*valX,     (k+1)*valY, l*valZ});
                            vert.push_back({h*valX,     k*valY,     (l+1)*valZ});
                            vert.push_back({(h+1)*valX, (k+1)*valY, l*valZ});
                            vert.push_back({h*valX,     (k+1)*valY, (l+1)*valZ});
                            vert.push_back({(h+1)*valX, k*valY,     (l+1)*valZ});
                        }
                    }
                }
            }else{
                // hk0
                for(int h=0; h < abs(hkl[0]); ++h){
                    for(int k=0; k < abs(hkl[1]); ++k){
                        vert.push_back({(h+1)*valX, k*valY, 0});
                        vert.push_back({(h+1)*valX, k*valY, 1});
                        vert.push_back({h*valX, (k+1)*valY, 0});
                        vert.push_back({(h+1)*valX, k*valY, 1});
                        vert.push_back({h*valX, (k+1)*valY, 0});
                        vert.push_back({h*valX, (k+1)*valY, 1});
                    }
                }
            }
        }else if(hasL){
            // h0l
            for(int h=0; h < abs(hkl[0]); ++h){
                for(int l=0; l < abs(hkl[2]); ++l){
                    vert.push_back({(h+1)*valX, 0, l*valZ});
                    vert.push_back({(h+1)*valX, 1, l*valZ});
                    vert.push_back({h*valX, 0, (l+1)*valZ});
                    vert.push_back({(h+1)*valX, 1, l*valZ});
                    vert.push_back({h*valX, 0, (l+1)*valZ});
                    vert.push_back({h*valX, 1, (l+1)*valZ});
                }
            }
        }else{
            // h00
            for(int h=1; h <= abs(hkl[0]); ++h){
                vert.push_back({h*valX, 0, 0});
                vert.push_back({h*valX, 1, 0});
                vert.push_back({h*valX, 0, 1});
                vert.push_back({h*valX, 1, 0});
                vert.push_back({h*valX, 0, 1});
                vert.push_back({h*valX, 1, 1});
            }
        }
    }else if(hasK){
        if(hasL){
            // 0lk
            for(int k=0; k < abs(hkl[1]); ++k){
                for(int l=0; l < abs(hkl[2]); ++l){
                    vert.push_back({0, (k+1)*valY, l*valZ});
                    vert.push_back({1, (k+1)*valY, l*valZ});
                    vert.push_back({0, k*valY, (l+1)*valZ});
                    vert.push_back({1, (k+1)*valY, l*valZ});
                    vert.push_back({0, k*valY, (l+1)*valZ});
                    vert.push_back({1, k*valY, (l+1)*valZ});
                }
            }
        }else{
            // 0l0
            for(int k=1; k <= abs(hkl[1]); ++k){
                vert.push_back({0, k*valY, 0});
                vert.push_back({1, k*valY, 0});
                vert.push_back({0, k*valY, 1});
                vert.push_back({1, k*valY, 0});
                vert.push_back({0, k*valY, 1});
                vert.push_back({1, k*valY, 1});
            }
        }
    }else if(hasL){
        // 00k
        for(int l=1; l <= abs(hkl[2]); ++l){
            vert.push_back({0, 0, l*valZ});
            vert.push_back({1, 0, l*valZ});
            vert.push_back({0, 1, l*valZ});
            vert.push_back({1, 0, l*valZ});
            vert.push_back({0, 1, l*valZ});
            vert.push_back({1, 1, l*valZ});
        }
    }
    for(size_t i=0; i<3; ++i){
        if(hkl[i] < 0){
            for(auto& v: vert){
                v[i] += 1;
            }
        }
    }
    auto cv = step->getCellVec() * step->getCellDim(CdmFmt::Bohr);
    for(auto& v: vert){
        v = v*cv;
    }
    return vert;
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
            active = curPlane->display;
            ui->pushButton->setChecked(active);
            if(active){
                master->addExtraData(&curPlane->gpu_data);
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
    }else if(change & (GuiChange::cell|GuiChange::settings)){
        if(curPlane){
            curPlane->gpu_data.update(settings.milCol.val);
            curPlane->gpu_data.update(mkVertices(curStep, curPlane->hkl));
        }
    }
}

void MillerWidget::updateIndex(int idx)
{
    if(curPlane){
        const auto* sender = QObject::sender();
        if(sender == ui->hSel){
            curPlane->hkl[0] = static_cast<int8_t>(idx);
        }else if(sender == ui->kSel){
            curPlane->hkl[1] = static_cast<int8_t>(idx);
        }else if(sender == ui->lSel){
            curPlane->hkl[2] = static_cast<int8_t>(idx);
        }else{
            throw Error("Unknown sender for HKL-plane index");
        }
        curPlane->gpu_data.update(mkVertices(curStep, curPlane->hkl));
    }
    if(active){
        triggerUpdate(GuiChange::extra);
    }
}

void MillerWidget::updateOffset(double off)
{
    if(curPlane){
        const auto* sender = QObject::sender();
        if(sender == ui->xOff){
            curPlane->offset[0] = static_cast<float>(off);
        }else if(sender == ui->yOff){
            curPlane->offset[1] = static_cast<float>(off);
        }else if(sender == ui->zOff){
            curPlane->offset[2] = static_cast<float>(off);
        }else{
            throw Error("Unknown sender for HKL-plane offset");
        }
        curPlane->gpu_data.update(curPlane->offset);
    }
    if(active){
        triggerUpdate(GuiChange::extra);
    }
}

void MillerWidget::on_pushButton_toggled(bool checked)
{
    active = checked;
    if(curPlane){
        curPlane->display = checked;
    }else if(active && (curPlane == nullptr)){
        auto hkl = std::array<int8_t,3>{
            static_cast<int8_t>(ui->hSel->value()),
            static_cast<int8_t>(ui->kSel->value()),
            static_cast<int8_t>(ui->lSel->value())};
        auto off = Vec{static_cast<float>(ui->xOff->value()),
                       static_cast<float>(ui->yOff->value()),
                       static_cast<float>(ui->zOff->value())};
        master->makeGLCurrent();
        auto tmp = planes.emplace(curStep, MillerPlane{
              true, hkl, off,
              GUI::MeshData{master->getGLGlobals(),
                            mkVertices(curStep, hkl),
                            off, settings.milCol.val, curStep}
              });
        curPlane = &tmp.first->second;
    }
    if(!curPlane){
        return;
    }
    if(active){
        master->addExtraData(&curPlane->gpu_data);
    }else{
        master->delExtraData(&curPlane->gpu_data);
    }
}
