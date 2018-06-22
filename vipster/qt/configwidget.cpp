#include "configwidget.h"
#include "ui_configwidget.h"
#include "iowrapper.h"

using namespace Vipster;

ConfigBase::ConfigBase(QWidget *parent):
    QWidget{parent}
{}

ConfigWidget::ConfigWidget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::ConfigWidget)
{
    ui->setupUi(this);
}

ConfigWidget::~ConfigWidget()
{
    delete ui;
}

void ConfigWidget::updateWidget(uint8_t change)
{
    if(change & Change::config){
        curConfig = master->curConfig;
        if(!curConfig){
            ui->configStack->setCurrentIndex(0);
            return;
        }
        if(auto pwConfig = dynamic_cast<IO::PWConfig*>(curConfig)){
            ui->configStack->setCurrentWidget(ui->PWWidget);
            ui->PWWidget->setConfig(pwConfig);
        }
        if(auto lmpConfig = dynamic_cast<IO::LmpConfig*>(curConfig)){
            ui->configStack->setCurrentWidget(ui->LmpWidget);
            ui->LmpWidget->setConfig(lmpConfig);
        }
    }
}

void ConfigWidget::registerConfig(IOFmt fmt, const std::string &name)
{
    ui->configSel->addItem(QString::fromStdString(
                               "(" + IOPlugins.at(fmt)->command +
                               ") " + name
                               ));
}

void ConfigWidget::on_configSel_currentIndexChanged(int index)
{
    master->setConfig(index);
}
