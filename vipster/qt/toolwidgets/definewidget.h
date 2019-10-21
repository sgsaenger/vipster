#ifndef DEFINEWIDGET_H
#define DEFINEWIDGET_H

#include <QWidget>
#include <QAction>

#include "molecule.h"
#include "../basewidget.h"
#include "../viewport.h"
#include "../common/seldata.h"

namespace Ui {
class DefineWidget;
}

class DefineWidget : public BaseWidget
{
    Q_OBJECT

public:
    explicit DefineWidget(QWidget *parent = nullptr);
    ~DefineWidget();
    void updateWidget(Vipster::GUI::change_t) override;

private slots:
    void on_newButton_clicked();
    void on_helpButton_clicked();
    void on_fromSelButton_clicked();

    void on_defTable_cellChanged(int row, int column);
    void on_defTable_itemSelectionChanged();

    void updateAction();
    void deleteAction();
    void toSelAction();

    void colButton_clicked();

private:
    void fillTable();

    Ui::DefineWidget *ui;
    Vipster::Step* curStep{nullptr};
    ViewPort::StepState *curState{nullptr};
    decltype(ViewPort::StepState::def)* defMap{nullptr};
    decltype(ViewPort::StepState::def)::iterator curSel;
    std::vector<std::shared_ptr<Vipster::GUI::SelData>> defList;
    QList<QAction*> contextActions;
};

#endif // DEFINEWIDGET_H
