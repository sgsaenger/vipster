#include "bondwidget.h"
#include "ui_bondwidget.h"
#include <QColorDialog>

using namespace Vipster;

BondColDelegate::BondColDelegate(QObject *parent)
    : QStyledItemDelegate{parent}
{}

QWidget *BondColDelegate::createEditor(QWidget *, const QStyleOptionViewItem &,
                                       const QModelIndex &index) const
{
    auto dialog = new QColorDialog{index.data(Qt::BackgroundRole).value<QColor>(), qApp->activeWindow()};
    dialog->setModal(Qt::WindowModality::ApplicationModal);
    return dialog;
}

void BondColDelegate::setModelData(QWidget *editor, QAbstractItemModel *model,
                                   const QModelIndex &index) const
{
    const auto& color = static_cast<QColorDialog*>(editor)->selectedColor();
    if(color.isValid()){
        model->setData(index, color, Qt::UserRole);
    }
}

bool BondColDelegate::eventFilter(QObject *object, QEvent *event)
{
    if(event->type() == QEvent::HideToParent){
        auto editor = static_cast<QColorDialog*>(object);
        emit commitData(editor);
        emit closeEditor(editor);
        return true;
    }
    return false;
}

void BondWidget::updateWidget(Vipster::GUI::change_t change)
{
    if (change & GUI::Change::atoms){
        bondModel.setStep(master->curStep);
    }
}

BondWidget::BondWidget(QWidget *parent) :
    BaseWidget(parent),
    ui(new Ui::BondWidget)
{
    ui->setupUi(this);
    ui->bondTable->setModel(&bondModel);
    ui->bondTable->setItemDelegateForColumn(3, new BondColDelegate{});
}

BondWidget::~BondWidget()
{
    delete ui;
}
