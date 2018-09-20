#include "lmpconfig.h"
#include "ui_lmpconfig.h"

using namespace Vipster;

LmpConfig::LmpConfig(QWidget *parent) :
    ConfigBase(parent),
    ui(new Ui::LmpConfig)
{
    ui->setupUi(this);
}

LmpConfig::~LmpConfig()
{
    delete ui;
}

void LmpConfig::setConfig(IO::BaseConfig *c)
{
    curConfig = dynamic_cast<IO::LmpConfig*>(c);
    if(!curConfig){
        throw Error("Invalid configuration preset");
    }
    QSignalBlocker block{this};
    ui->atomSel->setCurrentIndex(static_cast<int>(curConfig->style));
    ui->bondCheck->setCheckState(Qt::CheckState(curConfig->bonds*2));
    ui->angleCheck->setCheckState(Qt::CheckState(curConfig->angles*2));
    ui->dihedCheck->setCheckState(Qt::CheckState(curConfig->dihedrals*2));
    ui->impropCheck->setCheckState(Qt::CheckState(curConfig->impropers*2));
}

void LmpConfig::on_bondCheck_stateChanged(int arg1)
{
    curConfig->bonds = arg1/2;
}

void LmpConfig::on_angleCheck_stateChanged(int arg1)
{
    curConfig->angles = arg1/2;
}

void LmpConfig::on_dihedCheck_stateChanged(int arg1)
{
    curConfig->dihedrals = arg1/2;
}

void LmpConfig::on_impropCheck_stateChanged(int arg1)
{
    curConfig->impropers = arg1/2;
}

void LmpConfig::on_atomSel_currentIndexChanged(int index)
{
    curConfig->style = static_cast<IO::LmpConfig::AtomStyle>(index);
}
