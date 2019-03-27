#include "savefmtdialog.h"
#include "ui_savefmtdialog.h"
#include "mainwindow.h"

using namespace Vipster;

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
    auto* widget = ui->paramWidget;
    widget->clearParams();
    if(on){
        widget->setEnabled(true);
        const auto& mw = *static_cast<MainWindow*>(parentWidget());
        for(auto& p: mw.getParams()){
            if(p.first == fmt){
                widget->registerParam(p.second->copy());
            }
        }
        const auto& param_map = Vipster::params[fmt];
        for(const auto& p: param_map){
            widget->registerParam(p.second->copy());
        }
    }else{
        widget->setDisabled(true);
    }
}

void SaveFmtDialog::enableConfWidget(bool on)
{
    auto* widget = ui->configWidget;
    widget->clearConfigs();
    if(on){
        widget->setEnabled(true);
        const auto& mw = *static_cast<MainWindow*>(parentWidget());
        for(auto& p: mw.getConfigs()){
            if(p.first == fmt){
                widget->registerConfig(p.second->copy());
            }
        }
        const auto& conf_map = Vipster::configs[fmt];
        for(const auto& c: conf_map){
            widget->registerConfig(c.second->copy());
        }
    }else{
        widget->setDisabled(true);
    }
}

IO::BaseConfig* SaveFmtDialog::getConfig()
{
    return ui->configWidget->curConfig;
}

IO::BaseParam* SaveFmtDialog::getParam()
{
    return ui->paramWidget->curParam;
}
