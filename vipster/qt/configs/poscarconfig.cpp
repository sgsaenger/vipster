#include "poscarconfig.h"
#include "ui_poscarconfig.h"

using namespace Vipster;

PoscarConfig::PoscarConfig(QWidget *parent) :
    ConfigBase(parent),
    ui(new Ui::PoscarConfig)
{
    ui->setupUi(this);
}

PoscarConfig::~PoscarConfig()
{
    delete ui;
}

void PoscarConfig::setConfig(IO::BaseConfig *c)
{
    curConfig = dynamic_cast<IO::PoscarConfig*>(c);
    if(!curConfig){
        throw Error("Invalid configuration preset");
    }
    QSignalBlocker block{this};
    ui->fmtCombo->setCurrentIndex(curConfig->cartesian ? 1 : 0);
    ui->selCheck->setCheckState(curConfig->selective ?
                                Qt::CheckState::Checked : Qt::CheckState::Unchecked);
}

void PoscarConfig::on_fmtCombo_currentIndexChanged(int index)
{
    curConfig->cartesian = index;
}

void PoscarConfig::on_selCheck_toggled(bool checked)
{
    curConfig->selective = checked;
}
