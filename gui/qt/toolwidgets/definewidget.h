#ifndef DEFINEWIDGET_H
#define DEFINEWIDGET_H

#include <QWidget>
#include <QAction>
#include <type_traits>

#include "vipster/molecule.h"
#include "../vipsterapplication.h"
#include "seldata.h"

namespace Ui {
class DefineWidget;
}

class DefineWidget : public QWidget
{
    Q_OBJECT

public:
    explicit DefineWidget(QWidget *parent = nullptr);
    ~DefineWidget() override;

private slots:
    void createDefinition();
    void copySelToDefinition();

    void tableCellChanged(int row, int column);
    void tableSelectionChanged();

    void updateDefinition();
    void deleteDefinition();
    void copyDefToSelection();

    void changeColor();

private:
    void fillTable();

    Ui::DefineWidget *ui;
    decltype(Vipster::Application::StepState::definitions) *defMap{nullptr};
    std::remove_pointer_t<decltype(defMap)>::iterator curIt;
    Vipster::Step::const_selection &curSel();
    Vipster::SelectionFilter &curFilter();
    std::shared_ptr<Vipster::GUI::SelData> &curSelData();
    QList<QAction*> contextActions;
};

#endif // DEFINEWIDGET_H
