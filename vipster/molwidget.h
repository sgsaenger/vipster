#ifndef MOLWIDGET_H
#define MOLWIDGET_H

#include <QWidget>
#include <vector>
#include "mainwindow.h"
#include "molecule.h"

namespace Ui {
class MolWidget;
}

class MolWidget : public QWidget
{
    Q_OBJECT

public:
    explicit MolWidget(QWidget *parent = 0);
    ~MolWidget();
    void setStep(Vipster::Step *step);

signals:
    void stepChanged();
    void molChanged();

private slots:
    void on_cellDimBox_valueChanged(double cdm);
    void on_cellVecTable_cellChanged(int row, int column);
    void on_atomTable_cellChanged(int row, int column);

private:
    MainWindow *parent;
    Ui::MolWidget *ui;
    Vipster::Step *curStep{nullptr};
};

#endif // MOLWIDGET_H
