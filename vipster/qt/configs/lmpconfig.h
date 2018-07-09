#ifndef LMPCONFIG_H
#define LMPCONFIG_H

#include <QWidget>
#include "io/lmpinput/config.h"
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
    void setConfig(Vipster::IO::BaseConfig *c) override;

private slots:
    void on_bondCheck_stateChanged(int arg1);

    void on_angleCheck_stateChanged(int arg1);

    void on_dihedCheck_stateChanged(int arg1);

    void on_impropCheck_stateChanged(int arg1);

    void on_atomSel_currentIndexChanged(int index);

private:
    Ui::LmpConfig *ui;
    Vipster::IO::LmpConfig *curConfig{nullptr};
};

#endif // LMPCONFIG_H
