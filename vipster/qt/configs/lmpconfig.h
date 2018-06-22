#ifndef LMPCONFIG_H
#define LMPCONFIG_H

#include <QWidget>
#include "ioplugins/lmpinput.h"
#include "../configwidget.h"

namespace Ui {
class LmpConfig;
}

class LmpConfig : public ConfigBase
{
    Q_OBJECT

public:
    explicit LmpConfig(QWidget *parent = nullptr);
    ~LmpConfig() override;
    void setConfig(Vipster::BaseConfig *c) override;

private:
    Ui::LmpConfig *ui;
    Vipster::IO::LmpConfig *curConfig{nullptr};
};

#endif // LMPCONFIG_H
