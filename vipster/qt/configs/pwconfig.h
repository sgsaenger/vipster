#ifndef PWCONFIG_H
#define PWCONFIG_H

#include <QWidget>
#include "ioplugins/pwinput.h"
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

private:
    Ui::PWConfig *ui;
    Vipster::IO::PWConfig *curConfig{nullptr};
};

#endif // PWCONFIG_H
