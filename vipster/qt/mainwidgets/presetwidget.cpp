#include "presetwidget.h"
#include "ui_presetwidget.h"
#include "fileio.h"

#include <QMessageBox>
#include <QCheckBox>
#include <QComboBox>

using namespace Vipster;

PresetWidget::PresetWidget(QWidget *parent) :
    BaseWidget(parent),
    ui(new Ui::PresetWidget)
{
    ui->setupUi(this);
    ui->valueArea->hide();
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
                                  const IO::Preset &data)
{
    presets.emplace_back(name, data);
    ui->presetSel->addItem(QString::fromStdString(
                               "(" + data.getFmt()->command +
                               ") " + name
                               ));
    ui->presetSel->setCurrentIndex(ui->presetSel->count()-1);
}

void PresetWidget::on_presetSel_currentIndexChanged(int index)
{
    auto* layout = static_cast<QVBoxLayout*>(ui->valueArea->layout());
    // destroy old value-widgets
    while(auto *child = layout->takeAt(0)){
        while(auto *child_2 = child->layout()->takeAt(0)){
            layout->removeItem(child_2);
            child_2->widget()->deleteLater();
        }
        child->layout()->deleteLater();
    }
    // disable frame, show that nothing is loaded
    if(index<0){
        ui->noPreset->show();
        ui->valueArea->hide();
        curPreset = nullptr;
        return;
    }
    if(static_cast<size_t>(index) >= presets.size()){
        throw Error("Invalid IO preset selected");
    }
    // load values of new preset
    ui->noPreset->hide();
    ui->valueArea->show();
    auto& pair = presets.at(static_cast<size_t>(index));
    curPreset = &pair.second;
    for(auto& v: *curPreset){
        layout->insertLayout(0, new QHBoxLayout());
        auto *row = static_cast<QHBoxLayout*>(layout->children().back());
        row->addWidget(new QLabel{QString::fromStdString(v.first)+':'});
        switch(v.second.index()){
        case IO::Preset::i_bool:
        {
            auto box = new QCheckBox{};
            row->addWidget(box);
            box->setChecked(std::get<bool>(v.second));
            connect(box, &QCheckBox::stateChanged, [&v](int state){
                v.second = static_cast<bool>(state);
            });
            break;
        }
        case IO::Preset::i_enum:
        {
            auto box = new QComboBox{};
            auto &val = std::get<IO::Preset::i_enum>(v.second);
            row->addWidget(box);
            for(const auto& pair: val){
                box->addItem(QString::fromStdString(pair.second));
            }
            box->setCurrentIndex(val);
            connect(box, QOverload<int>::of(&QComboBox::currentIndexChanged), [&v](int i){
                std::get<IO::Preset::i_enum>(v.second) = i;
            });
            break;
        }
        default:
            row->addWidget(new QLabel{"(unsupported setting)"});
            break;
        }
    }
}

void PresetWidget::on_helpButton_clicked()
{
    QMessageBox::information(this, QString("About IO presets"), Vipster::IO::PresetsAbout);
}
