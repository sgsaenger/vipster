#include "xyzpreset.h"
#include "ui_xyzpreset.h"

using namespace Vipster;

XYZPreset::XYZPreset(QWidget *parent) :
    PresetBase(parent),
    ui(new Ui::XYZPreset)
{
    ui->setupUi(this);
}

XYZPreset::~XYZPreset()
{
    delete ui;
}

void XYZPreset::setPreset(IO::BasePreset *c)
{
    curPreset = dynamic_cast<IO::XYZPreset*>(c);
    if(!curPreset){
        throw Error("Invalid IO preset");
    }
    QSignalBlocker block{this};
    ui->dataSel->setCurrentIndex(static_cast<int>(curPreset->atomdata));
    ui->modeSel->setCurrentIndex(static_cast<int>(curPreset->filemode));
}

void XYZPreset::on_modeSel_currentIndexChanged(int index)
{
    curPreset->filemode = static_cast<IO::XYZPreset::Mode>(index);
}

void XYZPreset::on_dataSel_currentIndexChanged(int index)
{
    curPreset->atomdata = static_cast<IO::XYZPreset::Data>(index);
}
