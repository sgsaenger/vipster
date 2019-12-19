#include "pwpreset.h"
#include "ui_pwpreset.h"

using namespace Vipster;

PWPreset::PWPreset(QWidget *parent) :
    PresetBase(parent),
    ui(new Ui::PWPreset)
{
    ui->setupUi(this);
}

PWPreset::~PWPreset()
{
    delete ui;
}

void PWPreset::setPreset(IO::Preset *c)
{
    curPreset = c;
    if(curPreset->getFmt() != &IO::PWInput){
        throw Error("Invalid configuration preset");
    }
    ui->atomSel->setCurrentIndex(std::get<uint8_t>(curPreset->at("atoms")));
    ui->cellSel->setCurrentIndex(std::get<uint8_t>(curPreset->at("cell")));
}

void PWPreset::on_atomSel_currentIndexChanged(int index)
{
    curPreset->at("atoms") = static_cast<uint8_t>(index);
}

void PWPreset::on_cellSel_currentIndexChanged(int index)
{
    curPreset->at("cell") = static_cast<uint8_t>(index);
}
