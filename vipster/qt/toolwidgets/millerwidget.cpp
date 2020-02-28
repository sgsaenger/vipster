#include "../mainwindow.h"
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

MillerWidget::MillerPlane::MillerPlane(
        const Vipster::GUI::GlobalData& glob, std::vector<Face>&& faces,
        Vipster::Vec off, Vipster::Mat cell, Texture texture,
        const std::array<int8_t, 3> &hkl) :
    GUI::MeshData{glob, std::move(faces), off, cell, texture},
    hkl{hkl}
{}

std::vector<GUI::MeshData::Face> mkFaces(const std::array<int8_t,3>& hkl)
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
    return faces;
}

void MillerWidget::updateWidget(GUI::change_t change)
{
    if((change & GUI::stepChanged) == GUI::stepChanged){
        // set new pointers
        curStep = master->curStep;
        if(auto pos=planes.find(curStep); pos !=planes.end()){
            curPlane = pos->second;
        }else{
            curPlane = nullptr;
        }
        // set GUI state as needed
        if(curPlane){
            ui->hSel->setValue(curPlane->hkl[0]);
            ui->kSel->setValue(curPlane->hkl[1]);
            ui->lSel->setValue(curPlane->hkl[2]);
            ui->xOff->setValue(static_cast<double>(curPlane->offset[0]));
            ui->yOff->setValue(static_cast<double>(curPlane->offset[1]));
            ui->zOff->setValue(static_cast<double>(curPlane->offset[2]));
            // TODO: mask event?
            auto block = QSignalBlocker{ui->pushButton};
            ui->pushButton->setChecked(true);
        }else{
            ui->hSel->setValue(0);
            ui->kSel->setValue(0);
            ui->lSel->setValue(0);
            ui->xOff->setValue(0.);
            ui->yOff->setValue(0.);
            ui->zOff->setValue(0.);
            auto block = QSignalBlocker{ui->pushButton};
            ui->pushButton->setChecked(false);
        }
    }else if(curPlane){
        if(change & GUI::Change::settings){
            curPlane->update({{master->settings.milCol.val}, 1, 1});
        }
        if(change & GUI::Change::cell){
            curPlane->update(curStep->getCellVec()*curStep->getCellDim(CdmFmt::Bohr));
            curPlane->update(mkFaces(curPlane->hkl));
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
        curPlane->update(mkFaces(curPlane->hkl));
        triggerUpdate(GUI::Change::extra);
    }
}

void MillerWidget::updateOffset(double off)
{
    if(curPlane){
        const auto* sender = QObject::sender();
        if(sender == ui->xOff){
            curPlane->offset[0] = off;
        }else if(sender == ui->yOff){
            curPlane->offset[1] = off;
        }else if(sender == ui->zOff){
            curPlane->offset[2] = off;
        }else{
            throw Error("Unknown sender for HKL-plane offset");
        }
        curPlane->update(mkFaces(curPlane->hkl));
        triggerUpdate(GUI::Change::extra);
    }
}

void MillerWidget::on_pushButton_toggled(bool checked)
{
    if(checked){
        // create new Plane
        auto hkl = std::array<int8_t,3>{
            static_cast<int8_t>(ui->hSel->value()),
            static_cast<int8_t>(ui->kSel->value()),
            static_cast<int8_t>(ui->lSel->value())};
        auto off = Vec{ui->xOff->value(),
                       ui->yOff->value(),
                       ui->zOff->value()};
        planes.insert_or_assign(curStep,
            std::make_shared<MillerPlane>(
                master->globals, mkFaces(hkl), off,
                curStep->getCellVec()*curStep->getCellDim(CdmFmt::Bohr),
                MillerPlane::Texture{{master->settings.milCol.val}, 1, 1}, hkl
            ));
        curPlane = planes.at(curStep);
        master->curVP->addExtraData(curPlane, false);
    }else{
        // delete Plane
        master->curVP->delExtraData(curPlane, false);
    }
    triggerUpdate(GUI::Change::extra);
}
