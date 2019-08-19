#include "bondmodel.h"

using namespace Vipster;

BondModel::BondModel(BondWidget *parent)
    :parent{parent}
{}

void BondModel::setStep(Vipster::Step *curStep)
{
    beginResetModel();
    this->curStep = curStep;
    curBonds =  &curStep->getBonds(0, BondPolicy::Cell, BondFrequency::Never);
    endResetModel();
}

QVariant BondModel::headerData(int section, Qt::Orientation orientation, int role) const
{
    if(role != Qt::DisplayRole) return QVariant{};
    if(orientation == Qt::Horizontal){
        return colNames[section];
    }
    return section;
}

int BondModel::rowCount(const QModelIndex &parent) const
{
    if (parent.isValid() || !curBonds)
        return 0;

    return static_cast<int>(curBonds->size());
}

int BondModel::columnCount(const QModelIndex &parent) const
{
    if (parent.isValid() || !curBonds)
        return 0;

    return 3;
}

QVariant BondModel::data(const QModelIndex &index, int role) const
{
    if(!index.isValid() || !curBonds)
        return QVariant{};

    if(role == Qt::DisplayRole || role == Qt::EditRole){
        const auto& bond = curBonds->at(index.row());
        switch(index.column()){
        case 0:
            return QStringLiteral("%1-%2").arg(bond.at1).arg(bond.at2);
        case 1:
            return bond.dist;
        case 2:
            return "";
        }
    }
    return QVariant{};
}

Qt::ItemFlags BondModel::flags(const QModelIndex& index) const
{
    if (!index.isValid() || !curStep)
        return Qt::NoItemFlags;

    if(index.column() == 2)
        return Qt::ItemIsEnabled | Qt::ItemIsSelectable | Qt::ItemIsEditable;
    return Qt::ItemIsEnabled | Qt::ItemIsSelectable;
}
