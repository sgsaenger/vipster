#include "presetwidget.h"
#include "ui_presetwidget.h"
#include "io.h"

#include <QMessageBox>

using namespace Vipster;

PresetWidget::PresetWidget(QWidget *parent) :
    BaseWidget(parent),
    ui(new Ui::PresetWidget)
{
    ui->setupUi(this);
    formats = makePresetWidgets();
    for(auto& p: formats){
        ui->presetStack->addWidget(p.second);
    }
}

PresetWidget::~PresetWidget()
{
    delete ui;
}

void PresetWidget::clearPresets()
{
    presets.clear();
    ui->presetSel->clear();
}

void PresetWidget::registerPreset(const std::string& name,
                                  std::unique_ptr<Vipster::IO::BasePreset>&& data)
{
    presets.emplace_back(name, std::move(data));
    ui->presetSel->addItem(QString::fromStdString(
                               "(" + presets.back().second->getFmt()->command +
                               ") " + presets.back().first
                               ));
    ui->presetSel->setCurrentIndex(ui->presetSel->count()-1);
}

void PresetWidget::on_presetSel_currentIndexChanged(int index)
{
    if(index<0){
        ui->presetStack->setCurrentWidget(ui->NoCWidget);
        curPreset = nullptr;
        return;
    }
    if(static_cast<size_t>(index) >= presets.size()){
        throw Error("Invalid IO preset selected");
    }
    const auto& pair = presets.at(static_cast<size_t>(index));
    curFmt = pair.second->getFmt();
    auto pos = formats.find(curFmt);
    if(pos == formats.end()){
        throw Error("Invalid IO preset format");
    }
    curPreset = pair.second.get();
    ui->presetStack->setCurrentWidget(pos->second);
    pos->second->setPreset(curPreset);
}

void PresetWidget::on_helpButton_clicked()
{
    QMessageBox::information(this, QString("About IO presets"), Vipster::IO::PresetsAbout);
}
