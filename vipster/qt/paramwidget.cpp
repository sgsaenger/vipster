#include "paramwidget.h"
#include "ui_paramwidget.h"
#include "iowrapper.h"

using namespace Vipster;

ParamBase::ParamBase(QWidget *parent):
    QWidget(parent)
{}

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
        if(auto pwParam = dynamic_cast<IO::PWParam*>(curParam)){
            ui->paramStack->setCurrentIndex(1);
            ui->PWWidget->setParam(pwParam);
        }
    }
}

void ParamWidget::registerParam(IOFmt fmt, const std::string &name)
{
    ui->paramSel->addItem(QString::fromStdString(
                          "(" +  IOPlugins.at(fmt)->command +
                           ") " + name
                         ));
}

void ParamWidget::on_paramSel_currentIndexChanged(int index)
{
    master->setParam(index);
}
