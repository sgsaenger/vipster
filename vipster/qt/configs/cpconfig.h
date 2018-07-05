#ifndef CPCONFIG_H
#define CPCONFIG_H

#include <QWidget>
#include "io/cpmdinput/config.h"
#include "../configwidget.h"

namespace Ui {
class CPConfig;
}

class CPConfig : public ConfigBase
{
    Q_OBJECT

public:
    explicit CPConfig(QWidget *parent = nullptr);
    ~CPConfig() override;
    void setConfig(Vipster::BaseConfig *c) override;

private slots:
    void on_scaleSel_currentIndexChanged(int index);

    void on_angstromSel_stateChanged(int arg1);

private:
    Ui::CPConfig *ui;
    Vipster::IO::CPConfig *curConfig{nullptr};
};

#endif // CPCONFIG_H
