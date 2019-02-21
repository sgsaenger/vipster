#include "presetwidget.h"
#include "ui_presetwidget.h"
#include "io.h"

#include <QMessageBox>

using namespace Vipster;

PresetBase::PresetBase(QWidget *parent):
    QWidget{parent}
{}

PresetWidget::PresetWidget(QWidget *parent) :
    BaseWidget(parent),
    ui(new Ui::PresetWidget)
{
    ui->setupUi(this);
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
                               "(" + IOPlugins.at(presets.back().second->getFmt())->command +
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
    curPreset = pair.second.get();
    curFmt = curPreset->getFmt();
    switch (curFmt) {
    case IOFmt::XYZ:
        ui->presetStack->setCurrentWidget(ui->XYZWidget);
        ui->XYZWidget->setPreset(curPreset);
        break;
    case IOFmt::PWI:
        ui->presetStack->setCurrentWidget(ui->PWWidget);
        ui->PWWidget->setPreset(curPreset);
        break;
    case IOFmt::LMP:
        ui->presetStack->setCurrentWidget(ui->LmpWidget);
        ui->LmpWidget->setPreset(curPreset);
        break;
    case IOFmt::CPI:
        ui->presetStack->setCurrentWidget(ui->CPWidget);
        ui->CPWidget->setPreset(curPreset);
        break;
    default:
        throw Error("Invalid IO preset format");
    }
}

void PresetWidget::on_helpButton_clicked()
{
    QMessageBox::information(this, QString("About IO presets"), Vipster::IO::PresetsAbout);
}
