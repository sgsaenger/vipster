#include "cppreset.h"
#include "ui_cppreset.h"

using namespace Vipster;

CPPreset::CPPreset(QWidget *parent) :
    PresetBase(parent),
    ui(new Ui::CPPreset)
{
    ui->setupUi(this);
}

CPPreset::~CPPreset()
{
    delete ui;
}

void CPPreset::setPreset(IO::BasePreset *c)
{
    curPreset = c;
    if(curPreset->getFmt() != &IO::CPInput){
        throw Error("Invalid configuration preset");
    }
    QSignalBlocker block{this};
    ui->fmtSel->setCurrentIndex(std::get<uint>(curPreset->at("fmt")));
}

void CPPreset::on_fmtSel_currentIndexChanged(int index)
{
    curPreset->at("fmt") = static_cast<uint>(index);
}
