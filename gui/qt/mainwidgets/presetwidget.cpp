#include "presetwidget.h"
#include "ui_presetwidget.h"
#include "vipster/fileio.h"

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
                                  const Preset &data)
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
    curPreset = &presets.at(static_cast<size_t>(index)).second;
    for(auto& v: *curPreset){
        auto *row = new QHBoxLayout{};
        layout->addLayout(row);
        // setup label
        auto label = new QLabel{QString::fromStdString(v.first)+':'};
        row->addWidget(label);
        QWidget* widget{};
        // setup specific widget
        switch(v.second.first.index()){
        case Preset::i_bool:
        {
            auto box = new QCheckBox{};
            widget = box;
            box->setLayoutDirection(Qt::RightToLeft);
            box->setChecked(std::get<bool>(v.second.first));
            connect(box, &QCheckBox::stateChanged, [&v](int state){
                v.second.first = static_cast<bool>(state);
            });
            break;
        }
        case Preset::i_enum:
        {
            auto box = new QComboBox{};
            widget = box;
            auto &val = std::get<NamedEnum>(v.second.first);
            for(const auto& pair: val){
                box->addItem(QString::fromStdString(pair.second));
            }
            box->setCurrentIndex(val);
            connect(box, QOverload<int>::of(&QComboBox::currentIndexChanged), [&v](int i){
                std::get<NamedEnum>(v.second.first) = i;
            });
            break;
        }
        default:
            widget = new QLabel{"(unsupported setting)"};
            break;
        }
        // connect widget and setup tooltip
        row->addWidget(widget);
        label->setBuddy(widget);
        if(!v.second.second.empty()){
            label->setToolTip(v.second.second.c_str());
            widget->setToolTip(v.second.second.c_str());
        }
    }
}

void PresetWidget::on_helpButton_clicked()
{
    QMessageBox::information(this, QString("About IO presets"),
                             Vipster::PresetsAbout);
}
