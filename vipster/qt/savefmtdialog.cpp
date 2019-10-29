#include "savefmtdialog.h"
#include "ui_savefmtdialog.h"
#include "mainwindow.h"

using namespace Vipster;

SaveFmtDialog::SaveFmtDialog(const IO::Plugins &plugins, QWidget *parent) :
    QDialog(parent),
    ui(new Ui::SaveFmtDialog)
{
    ui->setupUi(this);
    for(const auto &plugin: plugins){
        if(plugin->writer){
            outFormats.push_back(plugin);
            ui->fmtSel->addItem(QString::fromStdString(plugin->name));
        }
    }
}

SaveFmtDialog::~SaveFmtDialog()
{
    delete ui;
}

void SaveFmtDialog::selFmt(int i)
{
    plugin = outFormats[static_cast<size_t>(i)];
    enableParamWidget(plugin->arguments & IO::Plugin::Args::Param);
    enableConfWidget(plugin->arguments & IO::Plugin::Args::Config);
}

void SaveFmtDialog::enableParamWidget(bool on)
{
    auto* widget = ui->paramWidget;
    widget->clearParams();
    if(on){
        widget->setEnabled(true);
        const auto& mw = *static_cast<MainWindow*>(parentWidget());
        for(const auto& p: mw.params.at(plugin)){
            widget->registerParam(p.second->copy());
        }
        for(auto& p: mw.getParams()){
            if(p.first == plugin){
                widget->registerParam(p.second->copy());
            }
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
        for(const auto& c: mw.configs.at(plugin)){
            widget->registerConfig(c.second->copy());
        }
        for(auto& p: mw.getConfigs()){
            if(p.first == plugin){
                widget->registerConfig(p.second->copy());
            }
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
