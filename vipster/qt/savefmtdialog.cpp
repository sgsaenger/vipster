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
    enablePresetWidget(iop->arguments & IO::Plugin::Args::Preset);
}

void SaveFmtDialog::enableParamWidget(bool on)
{
    auto* widget = ui->paramWidget;
    widget->clearParams();
    if(on){
        widget->setEnabled(true);
        const auto& mw = *static_cast<MainWindow*>(parentWidget());
        for(auto& p: mw.getParams()){
            if(p.second->getFmt() == fmt){
                widget->registerParam(p.first, p.second->copy());
            }
        }
        const auto& param_map = Vipster::params[fmt];
        for(const auto& p: param_map){
            widget->registerParam(p.first, p.second->copy());
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
        for(auto& p: mw.getPresets()){
            if(p.second->getFmt() == fmt){
                widget->registerPreset(p.first, p.second->copy());
            }
        }
        for(const auto& c: Vipster::presets[fmt]){
            widget->registerPreset(c.first, c.second->copy());
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
