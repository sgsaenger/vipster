#include "pwconfig.h"
#include "ui_pwconfig.h"

using namespace Vipster;

PWConfig::PWConfig(QWidget *parent) :
    ConfigBase(parent),
    ui(new Ui::PWConfig)
{
    ui->setupUi(this);
}

PWConfig::~PWConfig()
{
    delete ui;
}

void PWConfig::setConfig(IO::BaseConfig *c)
{
    curConfig = dynamic_cast<IO::PWConfig*>(c);
    if(!curConfig){
        throw Error("Invalid configuration preset");
    }
    ui->atomSel->setCurrentIndex(static_cast<int>(curConfig->atoms));
    ui->cellSel->setCurrentIndex(static_cast<int>(curConfig->cell));
}

void PWConfig::on_atomSel_currentIndexChanged(int index)
{
    curConfig->atoms = static_cast<IO::PWConfig::AtomFmt>(index);
}

void PWConfig::on_cellSel_currentIndexChanged(int index)
{
    curConfig->cell = static_cast<IO::PWConfig::CellFmt>(index);
}
