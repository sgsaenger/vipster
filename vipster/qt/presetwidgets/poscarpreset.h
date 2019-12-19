#ifndef POSCARCONFIG_H
#define POSCARCONFIG_H

#include <QWidget>
#include "io/plugins/poscar.h"
#include "../presetbase.h"

namespace Ui {
class PoscarPreset;
}

class PoscarPreset : public PresetBase
{
    Q_OBJECT

public:
    explicit PoscarPreset(QWidget *parent = nullptr);
    ~PoscarPreset();
    void setPreset(Vipster::IO::Preset *c) override;

private slots:
    void on_fmtCombo_currentIndexChanged(int index);
    void on_selCheck_toggled(bool checked);

private:
    Ui::PoscarPreset *ui;
    Vipster::IO::Preset *curPreset{nullptr};
};

#endif // POSCARCONFIG_H
