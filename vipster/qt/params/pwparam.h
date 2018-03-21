#ifndef PWIPARAM_H
#define PWIPARAM_H

#include <QWidget>
#include "ioplugins/pwinput.h"

namespace Ui {
class PWParam;
}

class PWParam : public QWidget
{
    Q_OBJECT

public:
    explicit PWParam(QWidget *parent = nullptr);
    ~PWParam();
    void setParam(Vipster::IO::PWParam *p);

private:
    Ui::PWParam *ui;
    Vipster::IO::PWParam *curParam{nullptr};
};

#endif // PWIPARAM_H
