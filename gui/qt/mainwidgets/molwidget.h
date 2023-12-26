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
