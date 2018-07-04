#ifndef PWCONFIG_H
#define PWCONFIG_H

#include <QWidget>
#include "io/pwinput/pwinput.h"
#include "../configwidget.h"

namespace Ui {
class PWConfig;
}

class PWConfig : public ConfigBase
{
    Q_OBJECT

public:
    explicit PWConfig(QWidget *parent = nullptr);
    ~PWConfig() override;
    void setConfig(Vipster::BaseConfig *c) override;

private slots:
    void on_atomSel_currentIndexChanged(int index);

    void on_cellSel_currentIndexChanged(int index);

private:
    Ui::PWConfig *ui;
    Vipster::IO::PWConfig *curConfig{nullptr};
};

#endif // PWCONFIG_H
