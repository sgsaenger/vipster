#include "molmodel.h"
#include "molwidget.h"

using namespace Vipster;

MolModel::MolModel(MolWidget *parent)
    : QAbstractTableModel(parent), parent{parent}
{
}

void MolModel::setStep(Vipster::Step *step)
{
    beginResetModel();
    curStep = step;
    endResetModel();
}

void MolModel::setColumns(int cols)
{
    beginResetModel();
    if(cols){
        colMap.clear();
        if(cols & 0x1){
            colMap.push_back(0);
        }
        if(cols & 0x2){
            colMap.push_back(1);
            colMap.push_back(2);
            colMap.push_back(3);
        }
        if(cols & 0x4){
            colMap.push_back(4);
        }
        if(cols & 0x8){
            colMap.push_back(5);
            colMap.push_back(6);
            colMap.push_back(7);
        }
        if(cols & 0x10){
            colMap.push_back(8);
        }
        if(cols & 0x20){
            colMap.push_back(9);
            colMap.push_back(10);
            colMap.push_back(11);
        }
    }
    endResetModel();
}

QVariant MolModel::headerData(int section, Qt::Orientation orientation, int role) const
{
    if(role != Qt::DisplayRole) return QVariant{};
    if(orientation == Qt::Horizontal){
        return colNames[colMap[section]];
    }else{
        return section+1;
    }
}

int MolModel::rowCount(const QModelIndex &parent) const
{
    if (parent.isValid() || !curStep)
        return 0;

    return curStep->getNat();
}

int MolModel::columnCount(const QModelIndex &parent) const
{
    if (parent.isValid() || !curStep)
        return 0;

    return colMap.size();
}

QVariant MolModel::data(const QModelIndex &index, int role) const
{
    if (!index.isValid() || !curStep)
        return QVariant{};

    if(role == Qt::DisplayRole || role == Qt::EditRole){
        int col = colMap[index.column()];
        const auto& atom = curStep->at(index.row());
        switch(col){
        case 0:
            return atom.name.c_str();
        case 1:
        case 2:
        case 3:
            return atom.coord[col-1];
        case 4:
            return atom.properties->charge;
        case 5:
        case 6:
        case 7:
            return atom.properties->forces[col-5];
        }
    }else if(role == Qt::CheckStateRole){
        int col = colMap[index.column()];
        const auto& atom = curStep->at(index.row());
        // return bool*2 so true becomes CheckState::Checked
        switch(col){
        case 8:
            return 2*atom.properties->flags[AtomFlag::Hidden];
        case 9:
            return 2*atom.properties->flags[AtomFlag::FixX];
        case 10:
            return 2*atom.properties->flags[AtomFlag::FixY];
        case 11:
            return 2*atom.properties->flags[AtomFlag::FixZ];
        }
    }
    return QVariant{};
}

bool MolModel::setData(const QModelIndex &index, const QVariant &value, int role)
{
    if (data(index, role) != value) {
        if(role == Qt::EditRole){
            int col = colMap[index.column()];
            auto atom = curStep->at(index.row());
            switch(col){
            case 0:
                atom.name = value.toString().toStdString();
                break;
            case 1:
            case 2:
            case 3:
                atom.coord[col-1] = value.toFloat();
                break;
            case 4:
                atom.properties->charge = value.toFloat();
                break;
            case 5:
            case 6:
            case 7:
                atom.properties->forces[col-5] = value.toFloat();
                break;
            }
        }else if(role == Qt::CheckStateRole){
            int col = colMap[index.column()];
            auto atom = curStep->at(index.row());
            switch(col){
            case 8:
                atom.properties->flags[AtomFlag::Hidden] = value==Qt::Checked;
                break;
            case 9:
                atom.properties->flags[AtomFlag::FixX] = value==Qt::Checked;
                break;
            case 10:
                atom.properties->flags[AtomFlag::FixY] = value==Qt::Checked;
                break;
            case 12:
                atom.properties->flags[AtomFlag::FixZ] = value==Qt::Checked;
                break;
            }
        }else{
            return false;
        }
        emit dataChanged(index, index, QVector<int>() << role);
        parent->triggerUpdate(GuiChange::atoms);
        return true;
    }
    return false;
}

Qt::ItemFlags MolModel::flags(const QModelIndex &index) const
{
    if (!index.isValid() || !curStep)
        return Qt::NoItemFlags;

    if(colMap[index.column()] >= 8)
        return Qt::ItemIsEnabled | Qt::ItemIsSelectable | Qt::ItemIsUserCheckable;
    return Qt::ItemIsEditable | Qt::ItemIsEnabled | Qt::ItemIsSelectable;
}
