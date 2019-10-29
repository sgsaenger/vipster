#include "configwidget.h"
#include "ui_configwidget.h"

#include <QMessageBox>

using namespace Vipster;

ConfigWidget::ConfigWidget(QWidget *parent) :
    BaseWidget(parent),
    ui(new Ui::ConfigWidget)
{
    ui->setupUi(this);
    formats = makeConfigWidgets();
    for(auto& p: formats){
        ui->configStack->addWidget(p.second);
    }
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
                               "(" + fmt->command +
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
    auto pos = formats.find(curFmt);
    if(pos == formats.end()){
        throw Error("Invalid parameter format");
    }
    curConfig = pair.second.get();
    ui->configStack->setCurrentWidget(pos->second);
    pos->second->setConfig(curConfig);
}

void ConfigWidget::on_helpButton_clicked()
{
    QMessageBox::information(this, QString("About IO-Configuration presets"), Vipster::IO::ConfigsAbout);
}
