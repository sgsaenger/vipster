#include <QCheckBox>
#include <QDoubleSpinBox>
#include <QColorDialog>
#include <QPushButton>
#include <QLineEdit>
#include <QComboBox>
#include <QSpacerItem>
#include "settingswidget.h"
#include "ui_settingswidget.h"
#include "vipsterapplication.h"

using namespace Vipster;

SettingsWidget::SettingsWidget(QWidget *parent) :
    BaseWidget(parent),
    ui(new Ui::SettingsWidget)
{
    ui->setupUi(this);
    registerSetting(&Settings::overlap);
    registerSetting(&Settings::atRadVdW);
    registerSetting(&Settings::atRadFac);
    registerSetting(&Settings::bondRad);
    registerSetting(&Settings::showCell);
    registerSetting(&Settings::antialias);
    registerSetting(&Settings::perspective);
    registerSetting(&Settings::rotCom);
    registerSetting(&Settings::animstep);
    registerSetting(&Settings::bgCol);
    registerSetting(&Settings::selCol);
    registerSetting(&Settings::milCol);
    registerSetting(&Settings::posCol);
    registerSetting(&Settings::negCol);
    ui->settingsLayout->addItem(new QSpacerItem(1,1, QSizePolicy::Preferred, QSizePolicy::Expanding), labels.size(), 0, 1, -1);
}

SettingsWidget::~SettingsWidget()
{
    delete ui;
}

template<>
QWidget* SettingsWidget::makeWidget(Setting<bool> Settings::* setting)
{
    auto* widget = new QCheckBox(ui->settingsContainer);
    widget->setCheckState(Qt::CheckState((vApp.config().settings.*setting).val * 2));
    connect(widget, &QCheckBox::toggled,
            [setting](bool checked){
                vApp.invokeOnConfig([](ConfigState &c, Setting<bool> Settings::* setting, bool checked){
                    (c.settings.*setting).val = checked;
                }, setting, checked);
            }
    );
    return widget;
}

template<>
QWidget* SettingsWidget::makeWidget(Setting<size_t> Settings::* setting)
{
    auto* widget = new QSpinBox(ui->settingsContainer);
    widget->setMaximum(10000);
    widget->setValue(static_cast<int>((vApp.config().settings.*setting).val));
    connect(widget, qOverload<int>(&QSpinBox::valueChanged),
            [setting](int newVal){
                vApp.invokeOnConfig([](ConfigState &c, Setting<size_t> Settings::* setting, size_t newVal){
                    (c.settings.*setting).val = newVal;
                }, setting, newVal);
            }
    );
    return widget;
}

template<>
QWidget* SettingsWidget::makeWidget(Setting<double> Settings::* setting)
{
    auto* widget = new QDoubleSpinBox(ui->settingsContainer);
    widget->setValue((vApp.config().settings.*setting).val);
    widget->setMinimum(0);
    connect(widget, qOverload<double>(&QDoubleSpinBox::valueChanged),
            [setting](double newVal){
                vApp.invokeOnConfig([](ConfigState &c, Setting<double> Settings::* setting, double newVal){
                    (c.settings.*setting).val = newVal;
                }, setting, newVal);
            }
    );
    return widget;
}

template<>
QWidget* SettingsWidget::makeWidget(Setting<std::string> Settings::* setting)
{
    auto* widget = new QLineEdit(ui->settingsContainer);
    widget->setText(QString::fromStdString((vApp.config().settings.*setting).val));
    connect(widget, &QLineEdit::editingFinished,
            [setting, widget](){
                vApp.invokeOnConfig([](ConfigState &c, Setting<std::string> Settings::* setting, const std::string& newVal){
                    (c.settings.*setting).val = newVal;
                }, setting, widget->text().toStdString());
            }
    );
    return widget;
}

template<>
QWidget* SettingsWidget::makeWidget(Setting<ColVec> Settings::* setting)
{
    auto* widget = new QPushButton(ui->settingsContainer);
    widget->setText("Select");
    auto origCol = (vApp.config().settings.*setting).val;
    widget->setStyleSheet(QString("background-color: rgb(%1,%2,%3)")
                          .arg(origCol[0]).arg(origCol[1]).arg(origCol[2]));
    connect(widget, &QPushButton::clicked,
            [setting, widget](){
                auto origCol = (vApp.config().settings.*setting).val;
                auto oldCol = QColor::fromRgb(origCol[0], origCol[1], origCol[2], origCol[3]);
                auto newCol = QColorDialog::getColor(oldCol, nullptr, QString{},
                                                     QColorDialog::ShowAlphaChannel);
                if(!newCol.isValid()){
                    return;
                }
                ColVec newVal = {static_cast<uint8_t>(newCol.red()),
                                 static_cast<uint8_t>(newCol.green()),
                                 static_cast<uint8_t>(newCol.blue()),
                                 static_cast<uint8_t>(newCol.alpha())};
                vApp.invokeOnConfig([](ConfigState &c, Setting<ColVec> Settings::* setting, const ColVec& newVal){
                    (c.settings.*setting).val = newVal;
                }, setting, newVal);
                widget->setStyleSheet(QString("background-color: %1").arg(newCol.name()));
            }
    );
    return widget;
}

template<typename T>
void SettingsWidget::registerSetting(Setting<T> Settings::* setting)
{
    const auto &settings = vApp.config().settings;

    labels.push_back(new QLabel(ui->settingsContainer));
    auto& label = labels.back();
    ui->settingsLayout->addWidget(label, static_cast<int>(labels.size())-1, 0);
    label->setText(QString::fromStdString((settings.*setting).name));

    widgets.push_back(makeWidget(setting));
    auto& widget = widgets.back();
    ui->settingsLayout->addWidget(widget, static_cast<int>(widgets.size())-1, 1, Qt::AlignRight);

    label->setBuddy(widget);
}
