#ifndef LMPCONFIG_H
#define LMPCONFIG_H

#include <QWidget>
#include "io/lmpinput/preset.h"
#include "../presetwidget.h"

namespace Ui {
class LmpPreset;
}

class LmpPreset : public PresetBase
{
    Q_OBJECT

public:
    explicit LmpPreset(QWidget *parent = nullptr);
    ~LmpPreset() override;
    void setPreset(Vipster::IO::BasePreset *c) override;

private slots:
    void on_bondCheck_stateChanged(int arg1);

    void on_angleCheck_stateChanged(int arg1);

    void on_dihedCheck_stateChanged(int arg1);

    void on_impropCheck_stateChanged(int arg1);

    void on_atomSel_currentIndexChanged(int index);

private:
    Ui::LmpPreset *ui;
    Vipster::IO::LmpPreset *curPreset{nullptr};
};

#endif // LMPCONFIG_H
