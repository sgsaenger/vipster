#include "bondwidget.h"
#include "ui_bondwidget.h"

using namespace Vipster;

void BondWidget::updateWidget(Vipster::guiChange_t change)
{
    if (change & GuiChange::atoms){
        bondModel.setStep(master->curStep);
    }
}

BondWidget::BondWidget(QWidget *parent) :
    BaseWidget(parent),
    ui(new Ui::BondWidget)
{
    ui->setupUi(this);
    ui->bondTable->setModel(&bondModel);
}

BondWidget::~BondWidget()
{
    delete ui;
}
