#include "bondmodel.h"
#include <QBrush>

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

    return 4;
}

QVariant BondModel::data(const QModelIndex &index, int role) const
{
    if(!index.isValid() || !curBonds)
        return QVariant{};

    if(role == Qt::DisplayRole || role == Qt::EditRole){
        const auto& bond = (*curBonds)[index.row()];
        switch(index.column()){
        case 0:
            return QStringLiteral("%1-%2").arg(bond.at1).arg(bond.at2);
        case 1:
            return bond.dist * invbohr;
        case 2:
            if (bond.type) {
                return bond.type->first.c_str();
            } else {
                return QStringLiteral("%1-%2").arg((*curStep)[bond.at1].name.c_str())
                                              .arg((*curStep)[bond.at2].name.c_str());
            }
        }
    }
    if(role == Qt::ForegroundRole && index.column() == 2){
        const auto& bond = (*curBonds)[index.row()];
        if(bond.type){
            return QBrush{QColor{0,0,0}};
        }else{
            return QBrush{QColor{128,128,128}};
        }
    }
    if(role == Qt::BackgroundRole && index.column() == 3){
        const auto& bond = (*curBonds)[index.row()];
        if(bond.type){
            const auto& col = bond.type->second;
            return QBrush{QColor{col[0],col[1],col[2]}};
        }else{
            return QBrush{QColor{255,255,255}};
        }
    }
    return QVariant{};
}

bool BondModel::setData(const QModelIndex &index, const QVariant &value, int role)
{
    if(role == Qt::EditRole){
        if(data(index, role) != value){
            if(index.column() == 2){
                curStep->setBondType(index.row(), value.toString().toStdString());
                return true;
            }
        }
    }else if(role == Qt::UserRole){
        const auto& color = value.value<QColor>();
        (*curBonds)[index.row()].type->second = {
                static_cast<uint8_t>(color.red()),
                static_cast<uint8_t>(color.green()),
                static_cast<uint8_t>(color.blue()),
                static_cast<uint8_t>(color.alpha())};
    }
    return false;
}

Qt::ItemFlags BondModel::flags(const QModelIndex& index) const
{
    if (!index.isValid() || !curStep)
        return Qt::NoItemFlags;

    if(index.column() == 2)
        return Qt::ItemIsEnabled | Qt::ItemIsEditable;
    if(index.column() == 3){
        if((*curBonds)[index.row()].type){
            return Qt::ItemIsEnabled | Qt::ItemIsEditable;
        }else{
            return Qt::NoItemFlags;
        }
    }
    return Qt::ItemIsEnabled;
}
