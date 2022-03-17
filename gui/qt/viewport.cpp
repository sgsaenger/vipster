#include "viewport.h"
#include "glwidget.h"
#include "mainwindow.h"
#include "ui_viewport.h"

using namespace Vipster;

ViewPort::ViewPort(MainWindow *parent, bool active) :
    QFrame(parent),
    ui(new Ui::ViewPort),
    master{parent},
    active{active}
{
    ui->setupUi(this);
    setFocusPolicy(Qt::StrongFocus);
    // try to create opengl-widget
    // TODO: catch error when no gl3.3 is available
    openGLWidget = new GLWidget{this, master->settings};
    ui->verticalLayout->insertWidget(1, openGLWidget, 1);
    setFocusProxy(openGLWidget);
    // connect timer for animation
    connect(&playTimer, &QTimer::timeout, ui->stepEdit, &QSpinBox::stepUp);
    // style buttons
    ui->firstStepButton->setIcon(style()->standardIcon(QStyle::SP_MediaSkipBackward));
    ui->preStepButton->setIcon(style()->standardIcon(QStyle::SP_MediaSeekBackward));
    ui->playButton->setIcon(style()->standardIcon(QStyle::SP_MediaPlay));
    ui->nextStepButton->setIcon(style()->standardIcon(QStyle::SP_MediaSeekForward));
    ui->lastStepButton->setIcon(style()->standardIcon(QStyle::SP_MediaSkipForward));
    ui->closeButton->setIcon(style()->standardIcon(QStyle::SP_TitleBarCloseButton));
    // fill mol-list
    for(const auto& mol: master->molecules){
        ui->molList->addItem(mol.name.c_str());
    }
}

ViewPort::ViewPort(const ViewPort &vp) :
    ViewPort{vp.master, false}
{
    //TODO: make this optional?
    moldata = vp.moldata;
    for(const auto& p: vp.stepdata){
        if(p.second.sel){
            if(stepdata[p.first].sel){
                *stepdata[p.first].sel = *p.second.sel;
            }else{
                stepdata[p.first].sel = std::make_unique<Step::selection>(*p.second.sel);
            }
        }
    }
    ui->molList->setCurrentIndex(vp.ui->molList->currentIndex());
    ui->stepEdit->setValue(vp.ui->stepEdit->value());
}

ViewPort::~ViewPort()
{
    delete ui;
}

void ViewPort::triggerUpdate(Vipster::GUI::change_t change)
{
    if(active || (curStep == master->curStep)){
        // if we are the active viewport or display the same step,
        // trigger global update
        master->updateWidgets(change);
    }else{
        // short-circuit if rest of GUI does not need to be updated
        // if necessary, make sure that bonds/overlaps are up to date
        if((change & GUI::Change::atoms) &&
           (master->stepdata[curStep].automatic_bonds ||
            master->settings.overlap.val)){
            curStep->generateBonds(!master->stepdata[curStep].automatic_bonds);
        }
        // trigger update in viewports that display the same step
        for(auto& vp: master->viewports){
            if(curStep == vp->curStep){
                vp->updateWidget(change);
            }
        }
    }
}

void ViewPort::updateWidget(GUI::change_t change)
{
    if(change & GUI::Change::cell){
        setMultEnabled(curStep->hasCell());
    }
    if(change & GUI::Change::extra){
        // remove extra-data that has been erased
        vpdata.extras.erase(
            std::remove_if(vpdata.extras.begin(), vpdata.extras.end(),
                           [](const auto &wp){return wp.expired();}),
            vpdata.extras.end());
        for(auto &sd: stepdata){
            sd.second.extras.erase(
                std::remove_if(sd.second.extras.begin(),
                               sd.second.extras.end(),
                               [](const auto &wp){return wp.expired();}),
                sd.second.extras.end());
        }
    }
    openGLWidget->updateWidget(change);
}

auto cmpExtraData(const std::shared_ptr<GUI::Data>& sp)
{
    return [&](const std::weak_ptr<GUI::Data>& wp){
        return (!wp.owner_before(sp)) && (!sp.owner_before(wp));
    };
}

void ViewPort::addExtraData(const std::shared_ptr<GUI::Data> &dat, bool global)
{
    auto &extras = global ? vpdata.extras : stepdata[curStep].extras;
    auto pos = std::find_if(extras.begin(), extras.end(), cmpExtraData(dat));
    if(pos == extras.end()){
        extras.push_back(dat);
    }
}

void ViewPort::delExtraData(const std::shared_ptr<GUI::Data> &dat, bool global)
{
    auto &extras = global ? vpdata.extras : stepdata[curStep].extras;
    auto pos = std::find_if(extras.begin(), extras.end(), cmpExtraData(dat));
    if(pos != extras.end()){
        extras.erase(pos);
    }
}

bool ViewPort::hasExtraData(const std::shared_ptr<GUI::Data> &dat, bool global)
{
    auto &extras = global ? vpdata.extras : stepdata[curStep].extras;
    return std::find_if(extras.begin(), extras.end(), cmpExtraData(dat)) != extras.end();
}

void ViewPort::registerMol(const std::string &name)
{
    ui->molList->addItem(name.c_str());
    if(active){
        // if we are the main viewport,
        // display new Molecule and make it active
        ui->molList->setCurrentIndex(ui->molList->count()-1);
    }
}

void ViewPort::setMol(int index)
{
    curMol = &*std::next(master->molecules.begin(), index);
    int nstep = static_cast<int>(curMol->getNstep());
    auto &curData = moldata[curMol];
    if(curData.curStep == 0){
        curData.curStep = nstep;
    }
    // set mol-selector
    QSignalBlocker selBlock{ui->molList};
    ui->molList->setCurrentIndex(index);
    // set mult manually
    QSignalBlocker xBlock{ui->xMultBox};
    ui->xMultBox->setValue(curData.mult[0]);
    QSignalBlocker yBlock{ui->yMultBox};
    ui->yMultBox->setValue(curData.mult[1]);
    QSignalBlocker zBlock{ui->zMultBox};
    ui->zMultBox->setValue(curData.mult[2]);
    openGLWidget->setMult(curData.mult);
    //Step-control
    ui->stepLabel->setText(QString::number(nstep));
    QSignalBlocker boxBlocker(ui->stepEdit);
    QSignalBlocker slideBlocker(ui->stepSlider);
    ui->stepEdit->setMaximum(nstep);
    ui->stepSlider->setMaximum(nstep);
    if(nstep == 1){
        ui->stepEdit->setDisabled(true);
        ui->stepSlider->setDisabled(true);
    }else{
        ui->stepEdit->setEnabled(true);
        ui->stepSlider->setEnabled(true);
    }
    // setStep will load the step and trigger the needed updates
    setStep(moldata[curMol].curStep, true);
}

void ViewPort::setStep(int i, bool setMol)
{
    curStep = &curMol->getStep(static_cast<size_t>(i-1));
    // ensure this step will be reused when mol is selected again
    moldata[curMol].curStep = i;
    // if step has not been accessed previously, set bonds to manual if bonds are already present
    if((master->stepdata.find(curStep) == master->stepdata.end()) &&
       !curStep->getBonds().empty()){
        master->stepdata[curStep].automatic_bonds = false;
    }
    // set widget's bond mode accordingly
    setBondMode(master->stepdata[curStep].automatic_bonds);
    // if no cell exists, disable mult-selectors
    setMultEnabled(curStep->hasCell());
    // if no previous selection exists in this viewport, create one, afterwards assign it
    auto& selRef = stepdata[curStep].sel;
    if(!selRef && curStep){
        selRef = std::make_unique<Step::selection>(curStep->select(SelectionFilter{}));
    }
    curSel = selRef.get();
    //Handle control-elements
    if(playTimer.isActive() && (i == static_cast<int>(curMol->getNstep()))){
        playTimer.stop();
        ui->playButton->setIcon(style()->standardIcon(QStyle::SP_MediaPlay));
    }
    QSignalBlocker boxBlocker(ui->stepEdit);
    QSignalBlocker slideBlocker(ui->stepSlider);
    ui->stepEdit->setValue(i);
    ui->stepSlider->setValue(i);
    if(i == 1){
        ui->preStepButton->setDisabled(true);
        ui->firstStepButton->setDisabled(true);
    }else{
        ui->preStepButton->setEnabled(true);
        ui->firstStepButton->setEnabled(true);
    }
    if(i == static_cast<int>(curMol->getNstep())){
        ui->nextStepButton->setDisabled(true);
        ui->lastStepButton->setDisabled(true);
    }else{
        ui->nextStepButton->setEnabled(true);
        ui->lastStepButton->setEnabled(true);
    }
    if(setMol){
        triggerUpdate(GUI::molChanged);
    }else{
        triggerUpdate(GUI::stepChanged);
    }
}

void ViewPort::setMult(int i)
{
    auto mult = moldata[curMol].mult;
    if(sender() == ui->xMultBox){ mult[0] = static_cast<uint8_t>(i); }
    else if(sender() == ui->yMultBox){ mult[1] = static_cast<uint8_t>(i); }
    else if(sender() == ui->zMultBox){ mult[2] = static_cast<uint8_t>(i); }
    // save mult if needed
    if (active) moldata[curMol].mult = mult;
    // trigger redraw
    openGLWidget->setMult(mult);
}

void ViewPort::setMultEnabled(bool b)
{
    ui->xMultBox->setEnabled(b);
    ui->yMultBox->setEnabled(b);
    ui->zMultBox->setEnabled(b);
}

void ViewPort::on_cameraGroup_buttonClicked(int i)
{
    openGLWidget->setCamera(i);
}

void ViewPort::setMouseMode(int i)
{
    if (i >= ui->mouseMode->count()) return;
    ui->mouseMode->setCurrentIndex(i);
}

void ViewPort::setBondMode(bool b)
{
    if(b){
        if(ui->mouseMode->currentIndex() == 3){
            ui->mouseMode->setCurrentIndex(0);
        }
        if(ui->mouseMode->count() == 4){
            ui->mouseMode->removeItem(3);
        }
    }else{
        if(ui->mouseMode->count() != 4){
            ui->mouseMode->addItem("Modify Bonds");
        }
    }
}

void ViewPort::stepBut(QAbstractButton *but)
{
    if(but == ui->firstStepButton){
        setStep(1);
    }else if(but == ui->lastStepButton){
        setStep(static_cast<int>(curMol->getNstep()));
    }else if(but == ui->playButton){
        if(playTimer.isActive()){
            playTimer.stop();
            ui->playButton->setIcon(style()->standardIcon(QStyle::SP_MediaPlay));
        }else{
            playTimer.start(static_cast<int>(master->settings.animstep.val));
            ui->playButton->setIcon(style()->standardIcon(QStyle::SP_MediaPause));
        }
    }
}

void ViewPort::on_mouseMode_currentIndexChanged(int i)
{
    openGLWidget->setMouseMode(static_cast<GLWidget::MouseMode>(i));
}

void ViewPort::on_closeButton_clicked()
{
    master->changeViewports(this, MainWindow::VP_CLOSE);
}

void ViewPort::on_vSplitButton_clicked()
{
    master->changeViewports(this, MainWindow::VP_VSPLIT);
}

void ViewPort::on_hSplitButton_clicked()
{
    master->changeViewports(this, MainWindow::VP_HSPLIT);
}

void ViewPort::makeActive(bool a)
{
    active = a;
    setFrameShadow(active ? QFrame::Sunken : QFrame::Raised);
}
