#ifndef SETTINGSWIDGET_H
#define SETTINGSWIDGET_H

#include <QWidget>
#include <QLabel>
#include "mainwindow.h"

namespace Ui {
class SettingsWidget;
}

class SettingsWidget : public QWidget, public BaseWidget
{
    Q_OBJECT

public:
    explicit SettingsWidget(QWidget *parent = nullptr);
    ~SettingsWidget();

private:
    template<typename T>
    void registerSetting(Vipster::Setting<T>& setting);
    template<typename T>
    QWidget* makeWidget(T& setting);
    std::vector<QLabel*> labels;
    std::vector<QWidget*> widgets;
    Ui::SettingsWidget *ui;
};

#endif // SETTINGSWIDGET_H
