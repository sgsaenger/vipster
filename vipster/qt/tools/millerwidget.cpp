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

std::vector<GUI::MeshData::Face> mkFaces(const std::array<int8_t,3>& hkl, Vec off)
{
    std::vector<GUI::MeshData::Face> faces{};
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
                            faces.push_back({{(h+1)*valX, k*valY,     l*valZ    },{},{}});
                            faces.push_back({{h*valX,     (k+1)*valY, l*valZ    },{},{}});
                            faces.push_back({{h*valX,     k*valY,     (l+1)*valZ},{},{}});
                            faces.push_back({{(h+1)*valX, (k+1)*valY, l*valZ    },{},{}});
                            faces.push_back({{h*valX,     (k+1)*valY, (l+1)*valZ},{},{}});
                            faces.push_back({{(h+1)*valX, k*valY,     (l+1)*valZ},{},{}});
                        }
                    }
                }
            }else{
                // hk0
                for(int h=0; h < abs(hkl[0]); ++h){
                    for(int k=0; k < abs(hkl[1]); ++k){
                        faces.push_back({{(h+1)*valX, k*valY, 0},{},{}});
                        faces.push_back({{(h+1)*valX, k*valY, 1},{},{}});
                        faces.push_back({{h*valX, (k+1)*valY, 0},{},{}});
                        faces.push_back({{(h+1)*valX, k*valY, 1},{},{}});
                        faces.push_back({{h*valX, (k+1)*valY, 0},{},{}});
                        faces.push_back({{h*valX, (k+1)*valY, 1},{},{}});
                    }
                }
            }
        }else if(hasL){
            // h0l
            for(int h=0; h < abs(hkl[0]); ++h){
                for(int l=0; l < abs(hkl[2]); ++l){
                    faces.push_back({{(h+1)*valX, 0, l*valZ},{},{}});
                    faces.push_back({{(h+1)*valX, 1, l*valZ},{},{}});
                    faces.push_back({{h*valX, 0, (l+1)*valZ},{},{}});
                    faces.push_back({{(h+1)*valX, 1, l*valZ},{},{}});
                    faces.push_back({{h*valX, 0, (l+1)*valZ},{},{}});
                    faces.push_back({{h*valX, 1, (l+1)*valZ},{},{}});
                }
            }
        }else{
            // h00
            for(int h=1; h <= abs(hkl[0]); ++h){
                faces.push_back({{h*valX, 0, 0},{},{}});
                faces.push_back({{h*valX, 1, 0},{},{}});
                faces.push_back({{h*valX, 0, 1},{},{}});
                faces.push_back({{h*valX, 1, 0},{},{}});
                faces.push_back({{h*valX, 0, 1},{},{}});
                faces.push_back({{h*valX, 1, 1},{},{}});
            }
        }
    }else if(hasK){
        if(hasL){
            // 0lk
            for(int k=0; k < abs(hkl[1]); ++k){
                for(int l=0; l < abs(hkl[2]); ++l){
                    faces.push_back({{0, (k+1)*valY, l*valZ},{},{}});
                    faces.push_back({{1, (k+1)*valY, l*valZ},{},{}});
                    faces.push_back({{0, k*valY, (l+1)*valZ},{},{}});
                    faces.push_back({{1, (k+1)*valY, l*valZ},{},{}});
                    faces.push_back({{0, k*valY, (l+1)*valZ},{},{}});
                    faces.push_back({{1, k*valY, (l+1)*valZ},{},{}});
                }
            }
        }else{
            // 0l0
            for(int k=1; k <= abs(hkl[1]); ++k){
                faces.push_back({{0, k*valY, 0},{},{}});
                faces.push_back({{1, k*valY, 0},{},{}});
                faces.push_back({{0, k*valY, 1},{},{}});
                faces.push_back({{1, k*valY, 0},{},{}});
                faces.push_back({{0, k*valY, 1},{},{}});
                faces.push_back({{1, k*valY, 1},{},{}});
            }
        }
    }else if(hasL){
        // 00k
        for(int l=1; l <= abs(hkl[2]); ++l){
            faces.push_back({{0, 0, l*valZ},{},{}});
            faces.push_back({{1, 0, l*valZ},{},{}});
            faces.push_back({{0, 1, l*valZ},{},{}});
            faces.push_back({{1, 0, l*valZ},{},{}});
            faces.push_back({{0, 1, l*valZ},{},{}});
            faces.push_back({{1, 1, l*valZ},{},{}});
        }
    }
    for(size_t i=0; i<3; ++i){
        if(hkl[i] < 0){
            for(auto& f: faces){
                f.pos[i] += 1;
            }
        }
    }
    for(auto& f: faces){
        f.pos += off;
    }
    return faces;
}

void MillerWidget::updateWidget(guiChange_t change)
{
    if((change & guiStepChanged) == guiStepChanged){
        //disable previous plane if necessary
        if(curPlane && curPlane->display){
            master->delExtraData(&curPlane->gpu_data);
        }
        // set new pointers
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
            if(curPlane->display){
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
        }
    }else if(curPlane){
        if(change & GuiChange::settings){
            curPlane->gpu_data.update({{settings.milCol.val}, 1, 1});
        }
        if(change & GuiChange::cell){
            curPlane->gpu_data.update(curStep->getCellVec()*curStep->getCellDim(CdmFmt::Bohr));
            curPlane->gpu_data.update(mkFaces(curPlane->hkl, curPlane->offset));
            if(curPlane->display != curStep->hasCell()){
                curPlane->display = curStep->hasCell();
                if(curPlane->display){
                    master->addExtraData(&curPlane->gpu_data);
                }else{
                    master->delExtraData(&curPlane->gpu_data);
                }
            }
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
        curPlane->gpu_data.update(mkFaces(curPlane->hkl, curPlane->offset));
        if(curPlane->display){
            triggerUpdate(GuiChange::extra);
        }
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
        curPlane->gpu_data.update(mkFaces(curPlane->hkl, curPlane->offset));
        if(curPlane->display){
            triggerUpdate(GuiChange::extra);
        }
    }
}

void MillerWidget::on_pushButton_toggled(bool checked)
{
    if(curPlane){
        curPlane->display = checked;
        if(checked){
            master->addExtraData(&curPlane->gpu_data);
        }else{
            master->delExtraData(&curPlane->gpu_data);
        }
    }else if(checked){
        auto hkl = std::array<int8_t,3>{
            static_cast<int8_t>(ui->hSel->value()),
            static_cast<int8_t>(ui->kSel->value()),
            static_cast<int8_t>(ui->lSel->value())};
        auto off = Vec{static_cast<float>(ui->xOff->value()),
                       static_cast<float>(ui->yOff->value()),
                       static_cast<float>(ui->zOff->value())};
        auto tmp = planes.emplace(curStep, MillerPlane{
              curStep->hasCell(), hkl, off,
              GUI::MeshData{master->getGLGlobals(),
                            mkFaces(hkl, off),
                            Vec{},
                            curStep->getCellVec()*curStep->getCellDim(CdmFmt::Bohr),
                            {{settings.milCol.val}, 1, 1}}
              });
        curPlane = &tmp.first->second;
        if(curPlane->display){
            master->addExtraData(&curPlane->gpu_data);
        }
    }
}
