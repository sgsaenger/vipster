#ifndef BONDWIDGET_H
#define BONDWIDGET_H

#include <QWidget>
#include <QStyledItemDelegate>
#include "../mainwindow.h"
#include "bondmodel.h"

namespace Ui {
class BondWidget;
}

class BondColDelegate : public QStyledItemDelegate
{
    Q_OBJECT
public:
    BondColDelegate(QObject *parent = nullptr);
    QWidget *createEditor(QWidget *parent, const QStyleOptionViewItem &option,
                          const QModelIndex &index) const override;
    void setModelData(QWidget *editor, QAbstractItemModel *model,
                      const QModelIndex &index) const override;
    bool eventFilter(QObject *object, QEvent *event) override;
};

class BondWidget : public BaseWidget
{
    Q_OBJECT

public:
    explicit BondWidget(QWidget *parent = nullptr);
    ~BondWidget() override;
    void updateWidget(Vipster::GUI::change_t change) override;

private:
    Ui::BondWidget *ui;
    BondModel bondModel{this};
};

#endif // BONDWIDGET_H
