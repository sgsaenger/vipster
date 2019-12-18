#ifndef XYZCONFIG_H
#define XYZCONFIG_H

#include <QWidget>
#include "io/xyz/plugin.h"
#include "../presetbase.h"

namespace Ui {
class XYZPreset;
}

class XYZPreset : public PresetBase
{
    Q_OBJECT

public:
    explicit XYZPreset(QWidget *parent = nullptr);
    ~XYZPreset();
    void setPreset(Vipster::IO::BasePreset *c) override;

private slots:
    void on_modeSel_currentIndexChanged(int index);

    void on_dataSel_currentIndexChanged(int index);

private:
    Ui::XYZPreset *ui;
    Vipster::IO::BasePreset *curPreset{nullptr};
};

#endif // XYZCONFIG_H
