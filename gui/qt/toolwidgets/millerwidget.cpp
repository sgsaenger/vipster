#include "../mainwindow.h"
#include "millerwidget.h"
#include "ui_millerwidget.h"
#include <QWindow>

using namespace Vipster;

static std::vector<GUI::MeshData::Face> mkFaces(const std::array<int8_t, 3>& hkl)
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


MillerWidget::MillerWidget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::MillerWidget)
{
    ui->setupUi(this);

    connect(&vApp, &MainWindow::activeStepChanged, this, [&](const Step &step){
        // set new pointers
        if(auto pos=planes.find(&vApp.curStep()); pos !=planes.end()){
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
    });

    connect(&vApp, &MainWindow::stepChanged, this, [&](const Step &step){
        if (curPlane && (&step == &vApp.curStep())) {
            curPlane->update(step.getCellVec() * step.getCellDim(AtomFmt::Bohr));
            curPlane->update(mkFaces(curPlane->hkl));
        }
    });

    connect(&vApp, &MainWindow::configChanged, this, [&](const ConfigState &cs){
        if (curPlane) {
            curPlane->update({{cs.settings.milCol.val}, 1, 1});
        }
    });

    connect(ui->pushButton, &QPushButton::clicked, this, [&](bool checked){
        if(checked){
            // create new Plane
            auto hkl = std::array<int8_t,3>{
                static_cast<int8_t>(ui->hSel->value()),
                static_cast<int8_t>(ui->kSel->value()),
                static_cast<int8_t>(ui->lSel->value())};
            auto off = Vec{ui->xOff->value(),
                           ui->yOff->value(),
                           ui->zOff->value()};
            planes.insert_or_assign(&vApp.curStep(),
                std::make_shared<MillerPlane>(
                    mkFaces(hkl), off,
                    vApp.curStep().getCellVec() * vApp.curStep().getCellDim(AtomFmt::Bohr),
                    MillerPlane::Texture{{vApp.config().settings.milCol.val}, 1, 1},
                    hkl
                ));
            curPlane = planes.at(&vApp.curStep());
            vApp.curVP->addExtraData(curPlane, false);
        }else{
            // delete Plane
            vApp.curVP->delExtraData(curPlane, false);
        }
    });

    connect(ui->hSel, &QSpinBox::valueChanged, this, &MillerWidget::updateIndex);
    connect(ui->kSel, &QSpinBox::valueChanged, this, &MillerWidget::updateIndex);
    connect(ui->lSel, &QSpinBox::valueChanged, this, &MillerWidget::updateIndex);
    connect(ui->xOff, &QDoubleSpinBox::valueChanged, this, &MillerWidget::updateOffset);
    connect(ui->yOff, &QDoubleSpinBox::valueChanged, this, &MillerWidget::updateOffset);
    connect(ui->zOff, &QDoubleSpinBox::valueChanged, this, &MillerWidget::updateOffset);
}

MillerWidget::~MillerWidget()
{
    delete ui;
}

MillerWidget::MillerPlane::MillerPlane(std::vector<Face>&& faces,
        Vipster::Vec off, Vipster::Mat cell, Texture texture,
        const std::array<int8_t, 3> &hkl) :
    GUI::MeshData{std::move(faces), off, cell, texture},
    hkl{hkl}
{}

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
    }
    vApp.curVP->updateState();
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
    }
    vApp.curVP->updateState();
}
