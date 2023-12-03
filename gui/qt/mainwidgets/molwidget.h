#ifndef MOLWIDGET_H
#define MOLWIDGET_H

#include <QWidget>
#include <QItemSelection>
#include "vipster/molecule.h"
#include "../basewidget.h"

namespace Ui {
class MolWidget;
}

class MolWidget : public BaseWidget
{
    Q_OBJECT

public:
    explicit MolWidget(QWidget *parent = nullptr);
    ~MolWidget() override;
    void updateWidget(Vipster::GUI::change_t change) override;

private slots:
    // kpoint slots
    void on_kFmtButton_clicked();
    void on_bands_stateChanged(int);
    void on_crystal_stateChanged(int);
    void mpg_change();
    void on_actionNew_K_Point_triggered();
    void on_actionDelete_K_Point_triggered();
    void on_discretetable_itemSelectionChanged();
    void on_discretetable_cellChanged(int row, int column);

    // pte slots
    void on_clearTableButton_clicked();
    void on_newElemButton_clicked();

private:
    void updateStep(Vipster::Step &step);

    void setActiveMol(Vipster::Molecule &mol);

    void fillKPoints(void);
    Ui::MolWidget *ui;
    int curKPoint{-1};
};

#endif // MOLWIDGET_H
