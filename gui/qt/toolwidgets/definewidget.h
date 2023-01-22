#ifndef DEFINEWIDGET_H
#define DEFINEWIDGET_H

#include <QWidget>
#include <QAction>
#include <type_traits>

#include "vipsterapplication.h"
#include "vipster/molecule.h"
#include "../basewidget.h"
#include "seldata.h"

namespace Ui {
class DefineWidget;
}

class DefineWidget : public BaseWidget
{
    Q_OBJECT

public:
    explicit DefineWidget(QWidget *parent = nullptr);
    ~DefineWidget() override;
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
    decltype(Vipster::Application::StepState::definitions) *defMap{nullptr};
    std::remove_pointer_t<decltype(defMap)>::iterator curIt;
    Vipster::Step::selection &curSel();
    Vipster::SelectionFilter &curFilter();
    std::shared_ptr<Vipster::GUI::SelData> &curSelData();
    QList<QAction*> contextActions;
};

#endif // DEFINEWIDGET_H
