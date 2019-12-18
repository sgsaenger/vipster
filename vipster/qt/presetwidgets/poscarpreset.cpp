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
    curPreset = c;
    if(curPreset->getFmt() != &IO::Poscar){
        throw Error("Invalid configuration preset");
    }
    QSignalBlocker block{this};
    ui->fmtCombo->setCurrentIndex(std::get<bool>(curPreset->at("cartesian")) ? 1 : 0);
    ui->selCheck->setChecked(std::get<bool>(curPreset->at("selective")));
}

void PoscarPreset::on_fmtCombo_currentIndexChanged(int index)
{
    curPreset->at("cartesian") = static_cast<uint>(index);
}

void PoscarPreset::on_selCheck_toggled(bool checked)
{
    curPreset->at("selective") = checked;
}
