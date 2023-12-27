#ifndef NEWELEMENT_H
#define NEWELEMENT_H

#include <QDialog>

namespace Ui {
class newelement;
}

class newelement : public QDialog
{
    Q_OBJECT

public:
    explicit newelement(bool isGlobal, QWidget *parent = nullptr);
    ~newelement();

private slots:
    void accept() override;

private:
    Ui::newelement *ui;
    bool isGlobal;
};

#endif // NEWELEMENT_H
