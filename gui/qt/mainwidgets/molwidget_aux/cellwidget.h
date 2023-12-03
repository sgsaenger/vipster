#ifndef CELLWIDGET_H
#define CELLWIDGET_H

#include <QWidget>
#include "vipster/step.h"

namespace Ui {
    class CellWidget;
}

class CellWidget : public QWidget
{
    Q_OBJECT

public:
    explicit CellWidget(QWidget *parent = nullptr);
    ~CellWidget();

private:
    void updateStep(Vipster::Step &step);
    void enableCell(bool enable);
    void cellFmtChanged(int idx);
    void cellDimChanged(double dim);
    void cellVecChanged(int row, int col);
    void applyToTrajectory();

    Ui::CellWidget *ui;
};

#endif // CELLWIDGET_H
