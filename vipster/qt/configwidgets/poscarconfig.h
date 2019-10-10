#ifndef POSCARCONFIG_H
#define POSCARCONFIG_H

#include <QWidget>
#include "io/poscar/config.h"
#include "../mainwidgets/configwidget.h"

namespace Ui {
class PoscarConfig;
}

class PoscarConfig : public ConfigBase
{
    Q_OBJECT

public:
    explicit PoscarConfig(QWidget *parent = nullptr);
    ~PoscarConfig();
    void setConfig(Vipster::IO::BaseConfig *c) override;

private slots:
    void on_fmtCombo_currentIndexChanged(int index);
    void on_selCheck_toggled(bool checked);

private:
    Ui::PoscarConfig *ui;
    Vipster::IO::PoscarConfig *curConfig{nullptr};
};

#endif // POSCARCONFIG_H
