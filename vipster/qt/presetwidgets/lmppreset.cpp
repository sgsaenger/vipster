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

void LmpPreset::setPreset(IO::BasePreset *c)
{
    curPreset = dynamic_cast<IO::LmpPreset*>(c);
    if(!curPreset){
        throw Error("Invalid configuration preset");
    }
    QSignalBlocker block{this};
    ui->atomSel->setCurrentIndex(static_cast<int>(curPreset->style));
    ui->bondCheck->setCheckState(Qt::CheckState(curPreset->bonds*2));
    ui->angleCheck->setCheckState(Qt::CheckState(curPreset->angles*2));
    ui->dihedCheck->setCheckState(Qt::CheckState(curPreset->dihedrals*2));
    ui->impropCheck->setCheckState(Qt::CheckState(curPreset->impropers*2));
}

void LmpPreset::on_bondCheck_stateChanged(int arg1)
{
    curPreset->bonds = arg1/2;
}

void LmpPreset::on_angleCheck_stateChanged(int arg1)
{
    curPreset->angles = arg1/2;
}

void LmpPreset::on_dihedCheck_stateChanged(int arg1)
{
    curPreset->dihedrals = arg1/2;
}

void LmpPreset::on_impropCheck_stateChanged(int arg1)
{
    curPreset->impropers = arg1/2;
}

void LmpPreset::on_atomSel_currentIndexChanged(int index)
{
    curPreset->style = static_cast<IO::LmpPreset::AtomStyle>(index);
}
