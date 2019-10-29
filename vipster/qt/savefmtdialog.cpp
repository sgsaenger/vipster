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
    enablePresetWidget(plugin->arguments & IO::Plugin::Args::Preset);
}

void SaveFmtDialog::enableParamWidget(bool on)
{
    auto* widget = ui->paramWidget;
    widget->clearParams();
    if(on){
        widget->setEnabled(true);
        const auto& mw = *static_cast<MainWindow*>(parentWidget());
        for(const auto& p: mw.params.at(plugin)){
            widget->registerParam(p.first, p.second->copy());
        }
        for(auto& p: mw.getParams()){
            if(p.second->getFmt() == plugin){
                widget->registerParam(p.first, p.second->copy());
            }
        }
    }else{
        widget->setDisabled(true);
    }
}

void SaveFmtDialog::enablePresetWidget(bool on)
{
    auto* widget = ui->presetWidget;
    widget->clearPresets();
    if(on){
        widget->setEnabled(true);
        const auto& mw = *static_cast<MainWindow*>(parentWidget());
        for(const auto& p: mw.presets.at(plugin)){
            widget->registerPreset(p.first, p.second->copy());
        }
        for(auto& p: mw.getPresets()){
            if(p.second->getFmt() == plugin){
                widget->registerPreset(p.first, p.second->copy());
            }
        }
    }else{
        widget->setDisabled(true);
    }
}

IO::BasePreset* SaveFmtDialog::getPreset()
{
    return ui->presetWidget->curPreset;
}

IO::BaseParam* SaveFmtDialog::getParam()
{
    return ui->paramWidget->curParam;
}
