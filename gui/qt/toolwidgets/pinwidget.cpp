#include "../mainwindow.h"
#include "pinwidget.h"
#include "ui_pinwidget.h"
#include <QMessageBox>

using namespace Vipster;

PinWidget::PinWidget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::PinWidget)
{
    ui->setupUi(this);

    connect(ui->addStep, &QPushButton::clicked, this, &PinWidget::addPin);
    connect(ui->delStep, &QPushButton::clicked, this, &PinWidget::delPin);
    connect(ui->stepList, &QListWidget::currentRowChanged, this, &PinWidget::selectPin);

    connect(ui->showCell, &QCheckBox::toggled, this, &PinWidget::toggleCell);
    connect(ui->showStep, &QCheckBox::toggled, this, &PinWidget::toggleVisible);
    connect(ui->repeatStep, &QCheckBox::toggled, this, &PinWidget::toggleRepeat);

    connect(ui->xOffset, &QDoubleSpinBox::valueChanged, this, &PinWidget::setOffset);
    connect(ui->yOffset, &QDoubleSpinBox::valueChanged, this, &PinWidget::setOffset);
    connect(ui->zOffset, &QDoubleSpinBox::valueChanged, this, &PinWidget::setOffset);
    connect(ui->xMultBox, &QSpinBox::valueChanged, this, &PinWidget::setMult);
    connect(ui->yMultBox, &QSpinBox::valueChanged, this, &PinWidget::setMult);
    connect(ui->zMultBox, &QSpinBox::valueChanged, this, &PinWidget::setMult);


    connect(ui->insertStep, &QPushButton::clicked, this, &PinWidget::insertStep);

    connect(ui->helpButton, &QPushButton::clicked, this, [&](){
        QMessageBox::information(this, QString("About pinning"),
            "Pinned Steps are drawn along the currently active Step.\n\n"
            "\"Repeating\" a Step means it is drawn in periodic repetitions "
            "of the active Step, i.e. with the periodicity of the active Step.\n"
            "Contrarily, specifying the multipliers for the pinned Step "
            "itself draws it with its own periodicity.\n"
            "Specifying the offset allows the pinned Step to be shifted against the active Step "
            "without having to modify its structure.\n\n"
            "Inserting a Step takes both offset and multipliers into account, "
            "and additionally performs cell vector fitting if requested.\n\n"
            "Cell vector fitting can be used e.g. to create commensurable super cells. "
            "If fitting is enabled for a direction, "
            "the lattice will be shrunken or stretched to "
            "match this dimension in the active Step."
        );
    });

    connect(&vApp, &MainWindow::activeStepChanged, [&](const Step& s){
        // only allow adding current step if it is not already pinned
        auto pos = std::find_if(pinnedSteps.begin(), pinnedSteps.end(), [&](const auto &v){
            return v->curStep == &s;
        });
        ui->addStep->setEnabled(pos == pinnedSteps.end());
        // update other gui state
        selectPin(ui->stepList->currentRow());
    });
    connect(&vApp, &MainWindow::configChanged, [&](const ConfigState &c){
        // update radii
        for (auto &dat: pinnedSteps) {
            dat->update(dat->curStep, c.settings.atRadVdW.val, c.settings.atRadFac.val, c.settings.bondRad.val);
        }
        if (!pinnedSteps.empty()) {
            vApp.curVP->updateState();
        }
    });
}

PinWidget::~PinWidget()
{
    delete ui;
}

PinWidget::PinnedStep::PinnedStep(const Step *step, const std::string& name,
                                  GUI::PBCVec mult)
    : GUI::StepData{step},
      name{name},
      mult{mult}
{}

void PinWidget::PinnedStep::draw(const Vec &off,
                                 const GUI::PBCVec &m,
                                 const Mat &cv, bool drawCell, void *context)
{
    Vec off_loc = off + this->offset;
    Mat cv_loc = curStep->getCellVec() * curStep->getCellDim(AtomFmt::Bohr);
    auto mult_loc = curStep->hasCell() ? mult : GUI::PBCVec{{1,1,1}};
    if(repeat){
        for(int x=0; x<m[0]; ++x){
            for(int y=0; y<m[1]; ++y){
                for(int z=0; z<m[2]; ++z){
                    StepData::draw(off_loc + x*cv[0] + y*cv[1] + z*cv[2],
                            mult_loc, cv_loc, drawCell & this->cell, context);
                }
            }
        }
    }else{
        StepData::draw(off_loc, mult_loc, cv_loc, drawCell & this->cell, context);
    }
}

void PinWidget::setMult(int i)
{
    if (!curPin) return;
    auto &mult = curPin->mult;
    if(sender() == ui->xMultBox){
        mult[0] = static_cast<uint8_t>(i);
    }else if(sender() == ui->yMultBox){
        mult[1] = static_cast<uint8_t>(i);
    }else if(sender() == ui->zMultBox){
        mult[2] = static_cast<uint8_t>(i);
    }
    vApp.curVP->updateState();
}

void PinWidget::setOffset(double d)
{
    if (!curPin) return;
    auto &off = curPin->offset;
    if(sender() == ui->xOffset){
        off[0] = d;
    }else if(sender() == ui->yOffset){
        off[1] = d;
    }else if(sender() == ui->zOffset){
        off[2] = d;
    }
    vApp.curVP->updateState();
}

void PinWidget::toggleCell(bool checked)
{
    if (!curPin) return;
    curPin->cell = checked;
    vApp.curVP->updateState();
}

void PinWidget::toggleRepeat(bool checked)
{
    if (!curPin) return;
    curPin->repeat = checked;
    vApp.curVP->updateState();
}

void PinWidget::delPin()
{
    // remove local infos
    ui->insertStep->setDisabled(true);
    if(curPin->curStep == &vApp.curStep()){
        ui->addStep->setEnabled(true);
    }
    auto pos2 = std::find(pinnedSteps.begin(), pinnedSteps.end(), curPin);
    vApp.curVP->delExtraData(*pos2, true);
    pinnedSteps.erase(pos2);
    delete ui->stepList->takeItem(ui->stepList->currentRow());
}

void PinWidget::addPin()
{
    ui->addStep->setDisabled(true);
    // add to list of steps
    pinnedSteps.push_back(std::make_shared<PinnedStep>(&vApp.curStep(),
        vApp.curMol().name + " (Step "
            + std::to_string(vApp.curVP->moldata[&vApp.curMol()].curStep) + ')',
        GUI::PBCVec{1,1,1}));
    const auto settings = vApp.config().settings;
    pinnedSteps.back()->update(pinnedSteps.back()->curStep,
                               settings.atRadVdW.val, settings.atRadFac.val,
                               settings.bondRad.val);
    ui->stepList->addItem(QString::fromStdString(pinnedSteps.back()->name));
    // enable in current viewport
    vApp.curVP->addExtraData(pinnedSteps.back(), true);
}

void PinWidget::selectPin(int currentRow)
{
    curPin = currentRow < 0 ? nullptr : pinnedSteps[currentRow];
    auto hasPin = static_cast<bool>(curPin);
    auto hasCell = hasPin ? curPin->curStep->hasCell() : false;
    ui->insertStep->setEnabled(hasPin ? curPin->curStep != &vApp.curStep() : false);
    ui->delStep->setEnabled(hasPin);
    ui->showStep->setEnabled(hasPin);
    ui->repeatStep->setEnabled(hasPin);
    ui->xOffset->setEnabled(hasPin);
    ui->yOffset->setEnabled(hasPin);
    ui->zOffset->setEnabled(hasPin);
    ui->showCell->setEnabled(hasCell);
    ui->xMultBox->setEnabled(hasCell);
    ui->yMultBox->setEnabled(hasCell);
    ui->zMultBox->setEnabled(hasCell);
    ui->xFit->setEnabled(hasCell);
    ui->yFit->setEnabled(hasCell);
    ui->zFit->setEnabled(hasCell);
    if(hasPin){
        QSignalBlocker block{ui->showStep};
        ui->showStep->setChecked(vApp.curVP->hasExtraData(curPin, true));
        QSignalBlocker block1{ui->showCell};
        ui->showCell->setChecked(curPin->cell);
        QSignalBlocker block2{ui->repeatStep};
        ui->repeatStep->setChecked(curPin->repeat);
        QSignalBlocker blockx{ui->xMultBox};
        ui->xMultBox->setValue(curPin->mult[0]);
        QSignalBlocker blocky{ui->yMultBox};
        ui->yMultBox->setValue(curPin->mult[1]);
        QSignalBlocker blockz{ui->zMultBox};
        ui->zMultBox->setValue(curPin->mult[2]);
        QSignalBlocker blockox{ui->xOffset};
        ui->xOffset->setValue(curPin->offset[0]);
        QSignalBlocker blockoy{ui->yOffset};
        ui->yOffset->setValue(curPin->offset[1]);
        QSignalBlocker blockoz{ui->zOffset};
        ui->zOffset->setValue(curPin->offset[2]);
    }
}

void PinWidget::insertStep()
{
    if (!curPin || (curPin->curStep == &vApp.curStep())) return;
    Step s = *curPin->curStep;
    s.asFmt(AtomFmt::Bohr).modShift(curPin->offset);
    std::array<bool,3> fit = {ui->xFit->isChecked(),
                              ui->yFit->isChecked(),
                              ui->zFit->isChecked()};
    if (s.hasCell() && (curPin->mult != GUI::PBCVec{{1,1,1}})) {
        const auto &m = curPin->mult;
        s.modMultiply(m[0], m[1], m[2]);
    }
    if (s.hasCell() && (fit != std::array<bool,3>{{false, false, false}})){
        auto fac = vApp.curStep().getCellDim(AtomFmt::Bohr) /
                s.getCellDim(AtomFmt::Bohr);
        auto cell = s.getCellVec();
        const auto& target = vApp.curStep().getCellVec();
        if (fit[0]) cell[0] = target[0] * fac;
        if (fit[1]) cell[1] = target[1] * fac;
        if (fit[2]) cell[2] = target[2] * fac;
        s.setCellVec(cell, true);
    }
    // immediately hide pinned step
    vApp.curVP->delExtraData(curPin, true);
    vApp.invokeOnStep(&Step::newAtoms<Step::atom_source>, s);
}

void PinWidget::toggleVisible(bool checked)
{
    if (!curPin) return;
    if (checked) {
        // insert into viewports extras
        vApp.curVP->addExtraData(curPin, true);
    }else{
        // remove from viewports extras
        vApp.curVP->delExtraData(curPin, true);
    }
}
