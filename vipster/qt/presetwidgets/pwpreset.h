#ifndef PWCONFIG_H
#define PWCONFIG_H

#include <QWidget>
#include "io/pwinput/plugin.h"
#include "../presetbase.h"

namespace Ui {
class PWPreset;
}

class PWPreset : public PresetBase
{
    Q_OBJECT

public:
    explicit PWPreset(QWidget *parent = nullptr);
    ~PWPreset() override;
    void setPreset(Vipster::IO::BasePreset *c) override;

private slots:
    void on_atomSel_currentIndexChanged(int index);

    void on_cellSel_currentIndexChanged(int index);

private:
    Ui::PWPreset *ui;
    Vipster::IO::PWPreset *curPreset{nullptr};
};

#endif // PWCONFIG_H
