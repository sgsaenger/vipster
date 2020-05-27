#ifndef NEWELEMENT_H
#define NEWELEMENT_H

#include <QDialog>
#include "vipster/periodictable.h"

namespace Ui {
class newelement;
}

class newelement : public QDialog
{
    Q_OBJECT

public:
    explicit newelement(Vipster::PeriodicTable& table, QWidget *parent = nullptr);
    ~newelement();

private slots:
    void accept() override;

private:
    Ui::newelement *ui;
    Vipster::PeriodicTable &table;
};

#endif // NEWELEMENT_H
