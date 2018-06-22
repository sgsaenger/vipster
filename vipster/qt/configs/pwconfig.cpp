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

void PWConfig::setConfig(BaseConfig *c)
{
    curConfig = static_cast<IO::PWConfig*>(c);
    ui->atomSel->setCurrentIndex(static_cast<int>(curConfig->atoms));
    ui->cellSel->setCurrentIndex(static_cast<int>(curConfig->cell));
}
