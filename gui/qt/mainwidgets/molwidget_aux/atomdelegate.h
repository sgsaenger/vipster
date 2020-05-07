#ifndef COORDDELEGATE_H
#define COORDDELEGATE_H

#include <QStyledItemDelegate>

class AtomDelegate: public QStyledItemDelegate
{
    Q_OBJECT
public:
    AtomDelegate(QObject *parent = nullptr);
    QWidget *createEditor(QWidget *parent, const QStyleOptionViewItem &option,
                          const QModelIndex &index) const override;
//    void setModelData(QWidget *editor, QAbstractItemModel *model,
//                      const QModelIndex &index) const override;
//    bool eventFilter(QObject *object, QEvent *event) override;
};

#endif // COORDDELEGATE_H
