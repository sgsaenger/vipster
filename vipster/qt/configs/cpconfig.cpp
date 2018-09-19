#include "cpconfig.h"
#include "ui_cpconfig.h"

using namespace Vipster;

CPConfig::CPConfig(QWidget *parent) :
    ConfigBase(parent),
    ui(new Ui::CPConfig)
{
    ui->setupUi(this);
}

CPConfig::~CPConfig()
{
    delete ui;
}

void CPConfig::setConfig(IO::BaseConfig *c)
{
    curConfig = dynamic_cast<IO::CPConfig*>(c);
    if(!curConfig){
        throw Error("Invalid configuration preset");
    }
    QSignalBlocker block{this};
    ui->angstromSel->setCheckState(Qt::CheckState(curConfig->angstrom*2));
    ui->scaleSel->setCurrentIndex(static_cast<int>(curConfig->scale));
}

void CPConfig::on_scaleSel_currentIndexChanged(int index)
{
    curConfig->scale = static_cast<IO::CPConfig::Scale>(index);
}

void CPConfig::on_angstromSel_stateChanged(int arg1)
{
    curConfig->angstrom = arg1/2;
}
