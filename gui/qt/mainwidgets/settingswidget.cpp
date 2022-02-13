#include <QCheckBox>
#include <QDoubleSpinBox>
#include <QColorDialog>
#include <QPushButton>
#include <QLineEdit>
#include <QComboBox>
#include <QSpacerItem>
#include "settingswidget.h"
#include "ui_settingswidget.h"
#include "../mainwindow.h"

using namespace Vipster;

SettingsWidget::SettingsWidget(QWidget *parent) :
    BaseWidget(parent),
    ui(new Ui::SettingsWidget)
{
    ui->setupUi(this);
    registerSetting(master->settings.overlap);
    registerSetting(master->settings.atRadVdW);
    registerSetting(master->settings.atRadFac);
    registerSetting(master->settings.bondRad);
    registerSetting(master->settings.showCell);
    registerSetting(master->settings.antialias);
    registerSetting(master->settings.perspective);
    registerSetting(master->settings.rotCom);
    registerSetting(master->settings.animstep);
    registerSetting(master->settings.bgCol);
    registerSetting(master->settings.selCol);
    registerSetting(master->settings.milCol);
    registerSetting(master->settings.posCol);
    registerSetting(master->settings.negCol);
    ui->settingsLayout->addItem(new QSpacerItem(1,1, QSizePolicy::Preferred, QSizePolicy::Expanding), labels.size(), 0, 1, -1);
}

SettingsWidget::~SettingsWidget()
{
    delete ui;
}

template<typename T>
QWidget* SettingsWidget::makeWidget(T&)
{
    return new QWidget(ui->settingsContainer);
}

template<>
QWidget* SettingsWidget::makeWidget(bool& setting)
{
    auto* widget = new QCheckBox(ui->settingsContainer);
    widget->setCheckState(Qt::CheckState(setting*2));
    connect(widget, &QCheckBox::toggled, this,
            [&setting, this](bool checked){
                setting = checked;
                triggerUpdate(GUI::Change::settings);
            }
    );
    return widget;
}

template<>
QWidget* SettingsWidget::makeWidget(size_t& setting)
{
    auto* widget = new QSpinBox(ui->settingsContainer);
    widget->setMaximum(10000);
    widget->setValue(static_cast<int>(setting));
    connect(widget, qOverload<int>(&QSpinBox::valueChanged), this,
            [&setting, this](int newVal){
                setting = static_cast<size_t>(newVal);
                triggerUpdate(GUI::Change::settings);
            }
    );
    return widget;
}

template<>
QWidget* SettingsWidget::makeWidget(double& setting)
{
    auto* widget = new QDoubleSpinBox(ui->settingsContainer);
    widget->setValue(setting);
    widget->setMinimum(0);
    connect(widget, qOverload<double>(&QDoubleSpinBox::valueChanged), this,
            [&setting, this](double newVal){
                setting = newVal;
                triggerUpdate(GUI::Change::settings);
            }
    );
    return widget;
}

template<>
QWidget* SettingsWidget::makeWidget(std::string& setting)
{
    auto* widget = new QLineEdit(ui->settingsContainer);
    widget->setText(QString::fromStdString(setting));
    connect(widget, &QLineEdit::editingFinished, this,
            [&setting, widget, this](){
                setting = widget->text().toStdString();
                triggerUpdate(GUI::Change::settings);
            }
    );
    return widget;
}

template<>
QWidget* SettingsWidget::makeWidget(Vipster::ColVec& setting)
{
    auto* widget = new QPushButton(ui->settingsContainer);
    widget->setText("Select");
    widget->setStyleSheet(QString("background-color: rgb(%1,%2,%3)")
                          .arg(setting[0]).arg(setting[1]).arg(setting[2]));
    connect(widget, &QPushButton::clicked, this,
            [&setting, widget, this](){
                auto oldCol = QColor::fromRgb(setting[0], setting[1], setting[2], setting[3]);
                auto newCol = QColorDialog::getColor(oldCol, this, QString{},
                                                     QColorDialog::ShowAlphaChannel);
                if(!newCol.isValid()){
                    return;
                }
                setting = {static_cast<uint8_t>(newCol.red()),
                           static_cast<uint8_t>(newCol.green()),
                           static_cast<uint8_t>(newCol.blue()),
                           static_cast<uint8_t>(newCol.alpha())};
                widget->setStyleSheet(QString("background-color: %1").arg(newCol.name()));
                triggerUpdate(GUI::Change::settings);
            }
    );
    return widget;
}

template<typename T>
void SettingsWidget::registerSetting(Setting<T>& setting)
{
    labels.push_back(new QLabel(ui->settingsContainer));
    auto& label = labels.back();
    ui->settingsLayout->addWidget(label, static_cast<int>(labels.size())-1, 0);
    label->setText(QString::fromStdString(setting.name));

    widgets.push_back(makeWidget(setting.val));
    auto& widget = widgets.back();
    ui->settingsLayout->addWidget(widget, static_cast<int>(widgets.size())-1, 1, Qt::AlignRight);

    label->setBuddy(widget);
}
