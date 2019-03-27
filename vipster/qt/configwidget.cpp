#include "configwidget.h"
#include "ui_configwidget.h"

#include <QMessageBox>

using namespace Vipster;

ConfigBase::ConfigBase(QWidget *parent):
    QWidget{parent}
{}

ConfigWidget::ConfigWidget(QWidget *parent) :
    BaseWidget(parent),
    ui(new Ui::ConfigWidget)
{
    ui->setupUi(this);
}

ConfigWidget::~ConfigWidget()
{
    delete ui;
}

void ConfigWidget::clearConfigs()
{
    configs.clear();
    ui->configSel->clear();
}

void ConfigWidget::registerConfig(std::unique_ptr<Vipster::IO::BaseConfig>&& data)
{
    auto fmt = data->getFmt();
    configs.emplace_back(fmt, std::move(data));
    ui->configSel->addItem(QString::fromStdString(
                               "(" + IOPlugins.at(fmt)->command +
                               ") " + configs.back().second->name
                               ));
    ui->configSel->setCurrentIndex(ui->configSel->count()-1);
}

void ConfigWidget::on_configSel_currentIndexChanged(int index)
{
    if(index<0){
        ui->configStack->setCurrentWidget(ui->NoCWidget);
        curConfig = nullptr;
        return;
    }
    if(static_cast<size_t>(index) >= configs.size()){
        throw Error("Invalid configuration preset selected");
    }
    const auto& pair = configs.at(static_cast<size_t>(index));
    curFmt = pair.first;
    curConfig = pair.second.get();
    switch (curFmt) {
    case IOFmt::XYZ:
        ui->configStack->setCurrentWidget(ui->XYZWidget);
        ui->XYZWidget->setConfig(curConfig);
        break;
    case IOFmt::PWI:
        ui->configStack->setCurrentWidget(ui->PWWidget);
        ui->PWWidget->setConfig(curConfig);
        break;
    case IOFmt::LMP:
        ui->configStack->setCurrentWidget(ui->LmpWidget);
        ui->LmpWidget->setConfig(curConfig);
        break;
    case IOFmt::CPI:
        ui->configStack->setCurrentWidget(ui->CPWidget);
        ui->CPWidget->setConfig(curConfig);
        break;
    default:
        throw Error("Invalid config format");
    }
}

void ConfigWidget::on_helpButton_clicked()
{
    QMessageBox::information(this, QString("About IO-Configuration presets"), Vipster::IO::ConfigsAbout);
}
