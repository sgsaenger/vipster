#include "doubledelegate.h"
#include <QDoubleSpinBox>

DoubleDelegate::DoubleDelegate(QObject *parent)
    : QStyledItemDelegate{parent}
{}

QWidget *DoubleDelegate::createEditor(QWidget *parent, const QStyleOptionViewItem &option,
                                     const QModelIndex &index) const
{
    if(index.data(Qt::EditRole).userType() == QVariant::Double){
        auto doubleSpinBox = new QDoubleSpinBox(parent);
        doubleSpinBox->setDecimals(8);
        doubleSpinBox->setRange(std::numeric_limits<double>::lowest(),
                               std::numeric_limits<double>::max());
        doubleSpinBox->setStepType(QDoubleSpinBox::StepType::AdaptiveDecimalStepType);
        return doubleSpinBox;
    }else{
        return QStyledItemDelegate::createEditor(parent, option, index);
    }
}
