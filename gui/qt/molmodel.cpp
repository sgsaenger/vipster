#include "molmodel.h"
#include "mainwindow.h"

MolModel::MolModel(QObject *parent)
    : QAbstractTableModel(parent)
{
}

QVariant MolModel::headerData(int section, Qt::Orientation orientation, int role) const
{
    if(role == Qt::DisplayRole){
        return QString{"Molecule"};
    }else{
        return {};
    }
}

//QModelIndex MolModel::index(int row, int column, const QModelIndex &parent) const
//{
//    // FIXME: Implement me!
//}

//QModelIndex MolModel::parent(const QModelIndex &index) const
//{
//    // FIXME: Implement me!
//}

int MolModel::rowCount(const QModelIndex &parent) const
{
    if (parent.isValid()) return 0;
    return vApp.molecules.size();
}

int MolModel::columnCount(const QModelIndex &parent) const
{
    if (!parent.isValid()) return 0;
    return 1;
}

QVariant MolModel::data(const QModelIndex &index, int role) const
{
    if (!index.isValid()) return {};
    if (role == Qt::DisplayRole) {
        auto it = vApp.molecules.begin();
        std::advance(it, index.row());
        return it->name.c_str();
    } else {
        return {};
    }
}
