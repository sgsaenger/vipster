#ifndef ATOMLIST_H
#define ATOMLIST_H

#include <QAction>
#include <QItemSelection>
#include <QWidget>
#include <optional>
#include "vipster/molecule.h"
#include "atommodel.h"

namespace Ui {
    class AtomList;
}

class AtomList : public QWidget
{
    Q_OBJECT

public:
    explicit AtomList(QWidget *parent = nullptr);
    ~AtomList();

private slots:
    void setActiveStep(const Vipster::Step &step, const Vipster::Step::selection &sel);
    void updateStep(const Vipster::Step &step);
    void updateSelection(const Vipster::Step::selection &sel);

private:
    void selectionChanged();
    void fmtButtonHandler();
    void fmtSelectionHandler(int index);

    Ui::AtomList *ui;

    AtomModel atomModel{};
    QList<QAction*> headerActions;
};

#endif // ATOMLIST_H
