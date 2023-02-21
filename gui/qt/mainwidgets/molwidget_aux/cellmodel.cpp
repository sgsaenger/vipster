#include "cellmodel.h"
#include "../molwidget.h"
#include <QMessageBox>

using namespace Vipster;

CellModel::CellModel(MolWidget *parent)
    :QAbstractTableModel(parent), parent{parent}
{}

void CellModel::setStep(Step *step)
{
    beginResetModel();
    curStep = step;
    vec = curStep->getCellVec();
    endResetModel();
}

Mat CellModel::getVec()
{
    return vec;
}

QVariant CellModel::headerData(int section, Qt::Orientation orientation, int role) const
{
    if(role != Qt::DisplayRole) return QVariant{};
    if(orientation == Qt::Horizontal){
        constexpr char colnames[] = {'x','y','z'};
        return QString{colnames[section]};
    }else{
        return section+1;
    }
}

int CellModel::rowCount(const QModelIndex &parent) const
{
    if(parent.isValid()) return 0;
    return 3;
}

int CellModel::columnCount(const QModelIndex &parent) const
{
    if(parent.isValid()) return 0;
    return 3;
}

QVariant CellModel::data(const QModelIndex &index, int role) const
{
    if(!index.isValid()) return {};
    if(role == Qt::DisplayRole || role == Qt::EditRole){
        return vec[index.row()][index.column()];
    }
    return {};
}

bool CellModel::setData(const QModelIndex &index, const QVariant &value, int role)
{
    if(data(index, role) != value){
        if(role == Qt::EditRole){
            vec[index.row()][index.column()] = value.toDouble();
            // early exit if not enabled
            if(!curStep->hasCell()) return true;
            // set vec of step
            auto scale = parent->scale();
            try{
                curStep->setCellVec(vec, scale);
            }catch(const Error &e){
                QMessageBox::critical(parent, "Error setting cell vectors", e.what());
                setStep(curStep);
                return false;
            }
            emit dataChanged(index, index, QVector<int>() << role);
            GUI::change_t change = GUI::Change::cell;
            if(scale) change |= GUI::Change::atoms;
            // short-circuit resetting the molModel on parent
            if(scale != (parent->ownStep->getFmt() == AtomFmt::Crystal)){
                parent->atomModel.setStep(parent->ownStep.get());
                // TODO
//                parent->setSelection();
            }
            parent->triggerUpdate(change);
            return true;
        }
    }
    return false;
}

Qt::ItemFlags CellModel::flags(const QModelIndex &index) const
{
    if(!index.isValid()) return Qt::NoItemFlags;
    return Qt::ItemIsEditable | Qt::ItemIsEnabled | Qt::ItemIsSelectable;
}
