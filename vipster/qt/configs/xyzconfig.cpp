#include "xyzconfig.h"
#include "ui_xyzconfig.h"

using namespace Vipster;

XYZConfig::XYZConfig(QWidget *parent) :
    ConfigBase(parent),
    ui(new Ui::XYZConfig)
{
    ui->setupUi(this);
}

XYZConfig::~XYZConfig()
{
    delete ui;
}

void XYZConfig::setConfig(IO::BaseConfig *c)
{
    curConfig = dynamic_cast<IO::XYZConfig*>(c);
    if(!curConfig){
        throw Error("Invalid configuration preset");
    }
    QSignalBlocker block{this};
    ui->dataSel->setCurrentIndex(static_cast<int>(curConfig->atomdata));
    ui->modeSel->setCurrentIndex(static_cast<int>(curConfig->filemode));
}

void XYZConfig::on_modeSel_currentIndexChanged(int index)
{
    curConfig->filemode = static_cast<IO::XYZConfig::Mode>(index);
}

void XYZConfig::on_dataSel_currentIndexChanged(int index)
{
    curConfig->atomdata = static_cast<IO::XYZConfig::Data>(index);
}
