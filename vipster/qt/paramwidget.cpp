#include "paramwidget.h"
#include "ui_paramwidget.h"
#include "iowrapper.h"

using namespace Vipster;

ParamWidget::ParamWidget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::ParamWidget)
{
    ui->setupUi(this);
}

ParamWidget::~ParamWidget()
{
    delete ui;
}

void ParamWidget::updateWidget(Vipster::Change change)
{
    if(change & Change::param){
        curParam = master->curParam;
        if(!curParam){
            ui->paramStack->setCurrentIndex(0);
            return;
        }
        auto pwParam = dynamic_cast<IO::PWParam*>(curParam);
        if(pwParam){
            ui->paramStack->setCurrentIndex(1);
            ui->PWWidget->setParam(pwParam);
        }
    }
}

void ParamWidget::registerParam(const std::string &name)
{
    ui->paramSel->addItem(name.c_str());
}

void ParamWidget::on_paramSel_currentIndexChanged(int index)
{
    master->setParam(index);
}
