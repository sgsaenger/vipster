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

class ConfigWidget : public QWidget
{
    Q_OBJECT

public:
    explicit ConfigWidget(QWidget *parent = nullptr);
    ~ConfigWidget() override;
    std::vector<std::pair< Vipster::IOFmt, std::unique_ptr<Vipster::BaseConfig>>> configs;
    void registerConfig(Vipster::IOFmt fmt,
                        std::unique_ptr<Vipster::BaseConfig>&& data);
    void clearConfigs();
    Vipster::BaseConfig* curConfig;

private slots:
    void on_configSel_currentIndexChanged(int index);

private:
    Ui::ConfigWidget *ui;
};

#endif // CONFIGWIDGET_H
