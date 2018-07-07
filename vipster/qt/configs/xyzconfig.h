#ifndef XYZCONFIG_H
#define XYZCONFIG_H

#include <QWidget>
#include "io/xyz/config.h"
#include "../configwidget.h"

namespace Ui {
class XYZConfig;
}

class XYZConfig : public ConfigBase
{
    Q_OBJECT

public:
    explicit XYZConfig(QWidget *parent = nullptr);
    ~XYZConfig();
    void setConfig(Vipster::IO::BaseConfig *c) override;

private slots:
    void on_modeSel_currentIndexChanged(int index);

    void on_dataSel_currentIndexChanged(int index);

private:
    Ui::XYZConfig *ui;
    Vipster::IO::XYZConfig *curConfig{nullptr};
};

#endif // XYZCONFIG_H
