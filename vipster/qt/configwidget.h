#ifndef CONFIGWIDGET_H
#define CONFIGWIDGET_H

#include <QWidget>
#include "mainwindow.h"

namespace Ui {
class ConfigWidget;
}

class ConfigBase: public QWidget
{
    Q_OBJECT

public:
    explicit ConfigBase(QWidget *parent = nullptr);
    virtual ~ConfigBase() = default;
    virtual void setConfig(Vipster::BaseConfig *c)=0;
};

class ConfigWidget : public QWidget, public BaseWidget
{
    Q_OBJECT

public:
    explicit ConfigWidget(QWidget *parent = nullptr);
    ~ConfigWidget() override;
    void updateWidget(uint8_t change) override;
    void registerConfig(Vipster::IOFmt, const std::string& name);

private slots:
    void on_configSel_currentIndexChanged(int index);

private:
    Ui::ConfigWidget *ui;
    Vipster::BaseConfig *curConfig;
};

#endif // CONFIGWIDGET_H
