#include "bonddelegate.h"
#include <QColorDialog>
#include <QApplication>

BondDelegate::BondDelegate(QObject *parent)
    : QStyledItemDelegate{parent}
{}

QWidget *BondDelegate::createEditor(QWidget *, const QStyleOptionViewItem &,
                                       const QModelIndex &index) const
{
    auto dialog = new QColorDialog{index.data(Qt::UserRole).value<QColor>(), qApp->activeWindow()};
    dialog->setWindowModality(Qt::WindowModality::ApplicationModal);
    dialog->setOption(QColorDialog::ShowAlphaChannel);
    return dialog;
}

void BondDelegate::setModelData(QWidget *editor, QAbstractItemModel *model,
                                   const QModelIndex &index) const
{
    const auto& color = static_cast<QColorDialog*>(editor)->selectedColor();
    if(color.isValid()){
        model->setData(index, color, Qt::UserRole);
    }
}

bool BondDelegate::eventFilter(QObject *object, QEvent *event)
{
    if(event->type() == QEvent::HideToParent){
        auto editor = static_cast<QColorDialog*>(object);
        emit commitData(editor);
        emit closeEditor(editor);
        return true;
    }
    return false;
}
