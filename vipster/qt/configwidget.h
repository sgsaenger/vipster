#ifndef CONFIGWIDGET_H
#define CONFIGWIDGET_H

#include <QWidget>
#include "basewidget.h"
#include "io.h"

namespace Ui {
class ConfigWidget;
}

class ConfigBase: public QWidget
{
    Q_OBJECT

public:
    explicit ConfigBase(QWidget *parent = nullptr);
    virtual ~ConfigBase() = default;
    virtual void setConfig(Vipster::IO::BaseConfig *c)=0;
};

class ConfigWidget : public BaseWidget
{
    Q_OBJECT

public:
    explicit ConfigWidget(QWidget *parent = nullptr);
    ~ConfigWidget() override;
    std::vector<std::pair< Vipster::IOFmt, std::unique_ptr<Vipster::IO::BaseConfig>>> configs;
    void registerConfig(std::unique_ptr<Vipster::IO::BaseConfig>&& data);
    void clearConfigs();
    Vipster::IOFmt curFmt;
    Vipster::IO::BaseConfig* curConfig{nullptr};

private slots:
    void on_configSel_currentIndexChanged(int index);

    void on_helpButton_clicked();

private:
    Ui::ConfigWidget *ui;
};

#endif // CONFIGWIDGET_H
