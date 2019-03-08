#ifndef DEFINEWIDGET_H
#define DEFINEWIDGET_H

#include <QWidget>
#include <QAction>
#include "../mainwindow.h"

namespace Ui {
class DefineWidget;
}

class DefineWidget : public BaseWidget
{
    Q_OBJECT

public:
    explicit DefineWidget(QWidget *parent = nullptr);
    ~DefineWidget();
    void updateWidget(Vipster::guiChange_t) override;

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
    struct GroupData{
        bool display;
        Vipster::ColVec color;
        Vipster::GUI::SelData gpu_data;
    };
    Vipster::Step* curStep{nullptr};
    std::map<std::string, Vipster::Step::selection>* defMap;
    std::map<Vipster::Step*, std::map<std::string, GroupData>> dataMap;
    int curSel{-1};
    std::vector<std::string> curNames;
    QList<QAction*> contextActions;
    constexpr static std::array<Vipster::ColVec, 5> colors{
        Vipster::ColVec{80, 0, 0, 80},
        Vipster::ColVec{0, 80, 0, 80},
        Vipster::ColVec{80, 80, 0, 80},
        Vipster::ColVec{80, 0, 80, 80},
        Vipster::ColVec{0, 80, 80, 80}
    };
};

#endif // DEFINEWIDGET_H
