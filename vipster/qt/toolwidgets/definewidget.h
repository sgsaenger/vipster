#ifndef DEFINEWIDGET_H
#define DEFINEWIDGET_H

#include <QWidget>
#include <QAction>
#include <type_traits>

#include "molecule.h"
#include "../basewidget.h"
#include "../viewport.h"
#include "../common/seldata.h"
#include "../mainwindow.h"

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
    decltype(MainWindow::StepState::definitions) *defMap{nullptr};
    std::remove_pointer_t<decltype(defMap)>::iterator curIt;
    // TODO: check if this aliasing works at all
    Vipster::Step::selection &curSel{std::get<0>(curIt->second)};
    Vipster::SelectionFilter &curFilter{std::get<1>(curIt->second)};
    std::shared_ptr<Vipster::GUI::SelData> &curSelData{std::get<2>(curIt->second)};
    QList<QAction*> contextActions;
};

#endif // DEFINEWIDGET_H
