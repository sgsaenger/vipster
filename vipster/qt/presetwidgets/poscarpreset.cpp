#include "poscarpreset.h"
#include "ui_poscarpreset.h"

using namespace Vipster;

PoscarPreset::PoscarPreset(QWidget *parent) :
    PresetBase(parent),
    ui(new Ui::PoscarPreset)
{
    ui->setupUi(this);
}

PoscarPreset::~PoscarPreset()
{
    delete ui;
}

void PoscarPreset::setPreset(IO::BasePreset *c)
{
    curPreset = dynamic_cast<IO::PoscarPreset*>(c);
    if(!curPreset){
        throw Error("Invalid configuration preset");
    }
    QSignalBlocker block{this};
    ui->fmtCombo->setCurrentIndex(curPreset->cartesian ? 1 : 0);
    ui->selCheck->setCheckState(curPreset->selective ?
                                Qt::CheckState::Checked : Qt::CheckState::Unchecked);
}

void PoscarPreset::on_fmtCombo_currentIndexChanged(int index)
{
    curPreset->cartesian = index;
}

void PoscarPreset::on_selCheck_toggled(bool checked)
{
    curPreset->selective = checked;
}
