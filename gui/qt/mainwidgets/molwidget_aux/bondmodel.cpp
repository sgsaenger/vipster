#include <QBrush>
#include "bondmodel.h"
#include "mainwindow.h"

using namespace Vipster;

void BondModel::reset()
{
    beginResetModel();
    curBonds = &vApp.curStep().getBonds();
    automatic_bonds = vApp.getState(vApp.curStep()).automatic_bonds;
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
    if (parent.isValid())
        return 0;

    if(automatic_bonds){
        return 3;
    }else{
        return 4;
    }
}

QVariant BondModel::data(const QModelIndex &index, int role) const
{
    if(!index.isValid() || !curBonds)
        return QVariant{};

    if(role == Qt::DisplayRole || role == Qt::EditRole){
        const auto& bond = (*curBonds)[index.row()];
        switch(index.column()){
        case 0:
            if(std::any_of(bond.diff.begin(), bond.diff.end(),
                           [](const auto&i)->bool{return i;})){
                return QStringLiteral("%1-%2 <%3,%4,%5>").arg(bond.at1).arg(bond.at2)
                        .arg(bond.diff[0]).arg(bond.diff[1]).arg(bond.diff[2]);
            }
            return QStringLiteral("%1-%2").arg(bond.at1).arg(bond.at2);
        case 1:
            return bond.dist * bohrrad;
        case 2:
            if (bond.type) {
                return bond.type->first.c_str();
            } else {
                const auto& curStep = vApp.curStep();
                const std::string& n1 = curStep[bond.at1].name;
                const std::string& n2 = curStep[bond.at2].name;
                return QStringLiteral("%1-%2").arg(std::min(n1, n2).c_str())
                                              .arg(std::max(n1, n2).c_str());
            }
        case 3:
            if(!bond.type){
                return "N/A";
            }
            break;
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
        }
    }
    if(role == Qt::UserRole && index.column() == 3){
        const auto& bond = (*curBonds)[index.row()];
        if(bond.type){
            const auto& col = bond.type->second;
            return QColor{col[0],col[1],col[2], col[3]};
        }else{
            return QColor{255,255,255, 255};
        }
    }
    return QVariant{};
}

bool BondModel::setData(const QModelIndex &index, const QVariant &value, int role)
{
    if(role == Qt::EditRole){
        if(data(index, role) != value){
            if(index.column() == 2){
                vApp.invokeOnStep(&Step::setBondType, index.row(), value.toString().toStdString());
                return true;
            }
        }
    }else if(role == Qt::UserRole){
        vApp.invokeOnStep([](Step &s, int i, const QColor &c){
            s.getBonds()[i].type->second = {
                static_cast<uint8_t>(c.red()),
                static_cast<uint8_t>(c.green()),
                static_cast<uint8_t>(c.blue()),
                static_cast<uint8_t>(c.alpha())};
        }, index.row(), value.value<QColor>());
        return true;
    }
    return false;
}

Qt::ItemFlags BondModel::flags(const QModelIndex& index) const
{
    if (!index.isValid())
        return Qt::NoItemFlags;

    if(index.column() == 2){
        if (!automatic_bonds) {
            return Qt::ItemIsEnabled | Qt::ItemIsEditable;
        } else {
            return Qt::ItemIsEnabled;
        }
    }
    if(index.column() == 3){
        if((*curBonds)[index.row()].type){
            return Qt::ItemIsEnabled | Qt::ItemIsEditable;
        }else{
            return Qt::NoItemFlags;
        }
    }
    return Qt::ItemIsEnabled;
}
