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
    ui->fmtSel->setCurrentIndex(static_cast<int>(curConfig->fmt));
}

void CPConfig::on_fmtSel_currentIndexChanged(int index)
{
    curConfig->fmt = static_cast<IO::CPConfig::AtomFmt>(index);
}
