#ifndef SETTINGSWIDGET_H
#define SETTINGSWIDGET_H

#include <QWidget>
#include <QLabel>
#include "vipster/settings.h"

namespace Ui {
class SettingsWidget;
}

class SettingsWidget : public QWidget
{
    Q_OBJECT

public:
    explicit SettingsWidget(QWidget *parent = nullptr);
    ~SettingsWidget();

private:
    template<typename T>
    void registerSetting(Vipster::Setting<T> Vipster::Settings::* setting);
    template<typename T>
    QWidget* makeWidget(Vipster::Setting<T> Vipster::Settings::* setting);
    std::vector<QLabel*> labels;
    std::vector<QWidget*> widgets;
    Ui::SettingsWidget *ui;
};

#endif // SETTINGSWIDGET_H
