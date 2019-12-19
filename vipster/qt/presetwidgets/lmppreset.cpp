#include "lmppreset.h"
#include "ui_lmppreset.h"

using namespace Vipster;

LmpPreset::LmpPreset(QWidget *parent) :
    PresetBase(parent),
    ui(new Ui::LmpPreset)
{
    ui->setupUi(this);
}

LmpPreset::~LmpPreset()
{
    delete ui;
}

void LmpPreset::setPreset(IO::Preset *c)
{
    curPreset = c;
    if(curPreset->getFmt() != &IO::LmpInput){
        throw Error("Invalid configuration preset");
    }
    QSignalBlocker block{this};
    ui->atomSel->setCurrentIndex(std::get<uint8_t>(curPreset->at("style")));
    ui->bondCheck->setChecked(std::get<bool>(curPreset->at("bonds")));
    ui->angleCheck->setChecked(std::get<bool>(curPreset->at("angles")));
    ui->dihedCheck->setChecked(std::get<bool>(curPreset->at("dihedrals")));
    ui->impropCheck->setChecked(std::get<bool>(curPreset->at("impropers")));
}

void LmpPreset::on_bondCheck_stateChanged(int arg1)
{
    curPreset->at("bonds") = (arg1/2) != 0;
}

void LmpPreset::on_angleCheck_stateChanged(int arg1)
{
    curPreset->at("angles") = (arg1/2) != 0;
}

void LmpPreset::on_dihedCheck_stateChanged(int arg1)
{
    curPreset->at("dihedrals") = (arg1/2) != 0;
}

void LmpPreset::on_impropCheck_stateChanged(int arg1)
{
    curPreset->at("impropers") = (arg1/2) != 0;
}

void LmpPreset::on_atomSel_currentIndexChanged(int index)
{
    curPreset->at("style") = static_cast<uint8_t>(index);
}
