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

void LmpConfig::setConfig(BaseConfig *c)
{
    curConfig = static_cast<IO::LmpConfig*>(c);
    ui->atomSel->setCurrentIndex(static_cast<int>(curConfig->style));
    ui->bondCheck->setCheckState(Qt::CheckState(curConfig->bonds*2));
    ui->angleCheck->setCheckState(Qt::CheckState(curConfig->angles*2));
    ui->dihedCheck->setCheckState(Qt::CheckState(curConfig->dihedrals*2));
    ui->impropCheck->setCheckState(Qt::CheckState(curConfig->impropers*2));
}
