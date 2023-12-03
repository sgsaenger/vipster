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
    void setActiveStep(Vipster::Step &step, Vipster::Step::selection &sel);
    void updateStep(Vipster::Step &step);
    void updateSelection(Vipster::Step::selection &sel);

private:
    void selectionChanged();
    void fmtButtonHandler();
    void fmtSelectionHandler(int index);

    Ui::AtomList *ui;

    std::optional<Vipster::Step::formatter> ownStep{};
    AtomModel atomModel{};
    QList<QAction*> headerActions;
};

#endif // ATOMLIST_H
