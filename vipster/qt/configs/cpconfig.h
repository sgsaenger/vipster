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
    void setConfig(Vipster::IO::BaseConfig *c) override;

private slots:
    void on_fmtSel_currentIndexChanged(int index);

private:
    Ui::CPConfig *ui;
    Vipster::IO::CPConfig *curConfig{nullptr};
};

#endif // CPCONFIG_H
