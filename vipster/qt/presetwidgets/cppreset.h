#ifndef CPCONFIG_H
#define CPCONFIG_H

#include <QWidget>
#include "io/cpmdinput/plugin.h"
#include "../presetbase.h"

namespace Ui {
class CPPreset;
}

class CPPreset : public PresetBase
{
    Q_OBJECT

public:
    explicit CPPreset(QWidget *parent = nullptr);
    ~CPPreset() override;
    void setPreset(Vipster::IO::BasePreset *c) override;

private slots:
    void on_fmtSel_currentIndexChanged(int index);

private:
    Ui::CPPreset *ui;
    Vipster::IO::BasePreset *curPreset{nullptr};
};

#endif // CPCONFIG_H
