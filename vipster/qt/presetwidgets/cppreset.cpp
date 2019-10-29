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
    curPreset = dynamic_cast<IO::CPPreset*>(c);
    if(!curPreset){
        throw Error("Invalid configuration preset");
    }
    QSignalBlocker block{this};
    ui->fmtSel->setCurrentIndex(static_cast<int>(curPreset->fmt));
}

void CPPreset::on_fmtSel_currentIndexChanged(int index)
{
    curPreset->fmt = static_cast<IO::CPPreset::AtomFmt>(index);
}
