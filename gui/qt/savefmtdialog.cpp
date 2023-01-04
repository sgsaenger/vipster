#include "savefmtdialog.h"
#include "ui_savefmtdialog.h"
#include "mainwindow.h"
#include "vipsterapplication.h"

using namespace Vipster;

SaveFmtDialog::SaveFmtDialog(const PluginList &plugins, QWidget *parent) :
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
    enableParamWidget(static_cast<bool>(plugin->makeParam));
    enablePresetWidget(static_cast<bool>(plugin->makePreset));
}

void SaveFmtDialog::enableParamWidget(bool on)
{
    auto* widget = ui->paramWidget;
    widget->clearParams();
    if(on){
        widget->setEnabled(true);
        const auto& mw = *static_cast<MainWindow*>(parentWidget());
        for(const auto& p: vApp.config.parameters.at(plugin)){
            widget->registerParam(p.first, p.second);
        }
        for(auto& p: mw.getParams()){
            if(p.second.getFmt() == plugin){
                widget->registerParam(p.first, p.second);
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
        for(const auto& p: vApp.config.presets.at(plugin)){
            widget->registerPreset(p.first, p.second);
        }
        for(auto& p: mw.getPresets()){
            if(p.second.getFmt() == plugin){
                widget->registerPreset(p.first, p.second);
            }
        }
    }else{
        widget->setDisabled(true);
    }
}

std::optional<Preset> SaveFmtDialog::getPreset()
{
    const auto& p = ui->presetWidget->curPreset;
    if(p){
        return *p;
    }else{
        return std::nullopt;
    }
}

std::optional<Parameter> SaveFmtDialog::getParam()
{
    const auto& p = ui->paramWidget->curParam;
    if(p){
        return *p;
    }else{
        return std::nullopt;
    }
}
