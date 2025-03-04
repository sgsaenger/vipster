#ifndef CELLMODWIDGET_H
#define CELLMODWIDGET_H

#include <QWidget>
#include "vipster/molecule.h"

namespace Ui {
class CellModWidget;
}

class CellModWidget : public QWidget
{
    Q_OBJECT

public:
    explicit CellModWidget(QWidget *parent = nullptr);
    ~CellModWidget() override;

private slots:
    void updateStep(const Vipster::Step &s);
    void on_wrapButton_clicked();
    void on_cropButton_clicked();
    void on_multButton_clicked();
    void on_alignButton_clicked();
    void on_reshapeButton_clicked();

private:
    Ui::CellModWidget *ui;
};

#endif // CELLMODWIDGET_H
