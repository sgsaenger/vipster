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

void PWPreset::setPreset(IO::BasePreset *c)
{
    curPreset = dynamic_cast<IO::PWPreset*>(c);
    if(!curPreset){
        throw Error("Invalid configuration preset");
    }
    ui->atomSel->setCurrentIndex(static_cast<int>(curPreset->atoms));
    ui->cellSel->setCurrentIndex(static_cast<int>(curPreset->cell));
}

void PWPreset::on_atomSel_currentIndexChanged(int index)
{
    curPreset->atoms = static_cast<IO::PWPreset::AtomFmt>(index);
}

void PWPreset::on_cellSel_currentIndexChanged(int index)
{
    curPreset->cell = static_cast<IO::PWPreset::CellFmt>(index);
}
