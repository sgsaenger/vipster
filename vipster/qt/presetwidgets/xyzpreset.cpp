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

void XYZPreset::setPreset(IO::Preset *c)
{
    curPreset = c;
    if(curPreset->getFmt() != &IO::XYZ){
        throw Error("Invalid configuration preset");
    }
    QSignalBlocker block{this};
    ui->dataSel->setCurrentIndex(std::get<uint8_t>(curPreset->at("atomdata")));
    ui->modeSel->setCurrentIndex(std::get<uint8_t>(curPreset->at("filemode")));
}

void XYZPreset::on_modeSel_currentIndexChanged(int index)
{
    curPreset->at("filemode") = static_cast<uint8_t>(index);
}

void XYZPreset::on_dataSel_currentIndexChanged(int index)
{
    curPreset->at("atomdata") = static_cast<uint8_t>(index);
}
