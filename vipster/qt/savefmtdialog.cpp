#include "savefmtdialog.h"
#include "ui_savefmtdialog.h"
#include "mainwindow.h"

using namespace Vipster;

constexpr int SaveFmtDialog::paramlist[];
constexpr int SaveFmtDialog::conflist[];

SaveFmtDialog::SaveFmtDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::SaveFmtDialog)
{
    ui->setupUi(this);
    for(const auto &iop: IOPlugins){
        if(iop.second->writer){
            outFormats.push_back(iop.first);
            ui->fmtSel->addItem(QString::fromStdString(iop.second->name));
        }
    }
}

SaveFmtDialog::~SaveFmtDialog()
{
    delete ui;
}

void SaveFmtDialog::selFmt(int i)
{
    fmt = outFormats[static_cast<size_t>(i)];
    const auto& iop = IOPlugins.at(fmt);
    enableParamWidget(iop->arguments & IO::Plugin::Args::Param);
    enableConfWidget(iop->arguments & IO::Plugin::Args::Config);
}

void SaveFmtDialog::enableParamWidget(bool on)
{
    ui->paramLabel->setDisabled(true);
    ui->paramSel->setDisabled(true);
    ui->paramWidget->setDisabled(true);
    if(on){
        ownParams.clear();
        auto vipRange = Vipster::params.equal_range(fmt);
        for(auto it=vipRange.first; it!=vipRange.second; ++it){
            ownParams.push_back(it->second->copy());
        }
        const auto& mw = *static_cast<MainWindow*>(parentWidget());
        for(auto& p: mw.params){
            ownParams.push_back(p.second->copy());
        }
    }
    if(on && ownParams.size()){
        ui->paramWidget->setCurrentIndex(paramlist[static_cast<int>(fmt)]);
        ui->paramSel->clear();
        for(auto& p: ownParams){
            ui->paramSel->addItem(QString::fromStdString(p->name));
        }
        ui->paramLabel->setEnabled(true);
        ui->paramSel->setEnabled(true);
        ui->paramWidget->setEnabled(true);
        on_paramSel_currentIndexChanged(0);
    }else{
        ui->paramSel->clear();
        ui->paramSel->addItem("None");
        ui->paramWidget->setCurrentIndex(0);
    }
}

void SaveFmtDialog::enableConfWidget(bool on)
{
    ui->configLabel->setDisabled(true);
    ui->configSel->setDisabled(true);
    ui->configWidget->setDisabled(true);
    if(on){
        ownConfigs.clear();
        // TODO: No default configs yet
//        auto vipRange = Vipster::config.equal_range(fmt);
//        for(auto it=vipRange.first; it!=vipRange.second; ++it){
//            ownParams.push_back(it->second->copy());
//        }
        const auto& mw = *static_cast<MainWindow*>(parentWidget());
        for(auto& p: mw.configs){
            ownConfigs.push_back(p.second->copy());
        }
    }
    if(on && ownConfigs.size()){
        ui->configWidget->setCurrentIndex(paramlist[static_cast<int>(fmt)]);
        ui->configSel->clear();
        for(auto& p: ownConfigs){
            ui->configSel->addItem(QString::fromStdString(p->name));
        }
        ui->configLabel->setEnabled(true);
        ui->configSel->setEnabled(true);
        ui->configWidget->setEnabled(true);
        on_configSel_currentIndexChanged(0);
    }else{
        ui->configSel->clear();
        ui->configSel->addItem("None");
        ui->configWidget->setCurrentIndex(0);
    }
}

void SaveFmtDialog::on_paramSel_currentIndexChanged(int index)
{
    if(ui->paramWidget->isEnabled()){
        param = ownParams[static_cast<size_t>(index)].get();
        static_cast<ParamBase*>(ui->paramWidget->
                                currentWidget())->
                setParam(param);
    }
}

void SaveFmtDialog::on_configSel_currentIndexChanged(int index)
{
    if(ui->configWidget->isEnabled()){
        config = ownConfigs[static_cast<size_t>(index)].get();
        //TODO: no config widget yet
//        static_cast<ConfigBase*>(ui->configWidget->
//                                 currentWidget())->
//                setConfig(config);
    }
}
