#include "atommodel.h"
#include "vipsterapplication.h"

using namespace Vipster;

enum colIds{
    c_name,
    c_x, c_y, c_z,
    c_charge,
    c_fx, c_fy, c_fz,
    c_hidden,
    c_fixx, c_fixy, c_fixz
};

void AtomModel::setStep(Step::formatter &&step)
{
    beginResetModel();
    curStep = step;
    endResetModel();
}

void AtomModel::update()
{
    endResetModel();
}

void AtomModel::setColumns(int cols)
{
    beginResetModel();
    if(cols){
        colMap.clear();
        if(cols & 0x10){
            colMap.push_back(c_hidden);
        }
        if(cols & 0x1){
            colMap.push_back(c_name);
        }
        if(cols & 0x2){
            colMap.push_back(c_x);
            colMap.push_back(c_y);
            colMap.push_back(c_z);
        }
        if(cols & 0x4){
            colMap.push_back(c_charge);
        }
        if(cols & 0x8){
            colMap.push_back(c_fx);
            colMap.push_back(c_fy);
            colMap.push_back(c_fz);
        }
        if(cols & 0x20){
            colMap.push_back(c_fixx);
            colMap.push_back(c_fixy);
            colMap.push_back(c_fixz);
        }
    }
    endResetModel();
}

QVariant AtomModel::headerData(int section, Qt::Orientation orientation, int role) const
{
    if(role != Qt::DisplayRole) return QVariant{};
    if(orientation == Qt::Horizontal){
        return colNames[colMap[section]];
    }else{
        return section;
    }
}

int AtomModel::rowCount(const QModelIndex &parent) const
{
    if (parent.isValid() || !curStep)
        return 0;

    return curStep->getNat();
}

int AtomModel::columnCount(const QModelIndex &parent) const
{
    if (parent.isValid() || !curStep)
        return 0;

    return colMap.size();
}

QVariant AtomModel::data(const QModelIndex &index, int role) const
{
    if (!index.isValid() || !curStep)
        return QVariant{};

    if(role == Qt::DisplayRole || role == Qt::EditRole){
        int col = colMap[index.column()];
        const auto& atom = curStep->at(index.row());
        switch(col){
        case c_name:
            return atom.name.c_str();
        case c_x:
        case c_y:
        case c_z:
            return atom.coord[col-1];
        case c_charge:
            return atom.properties->charge;
        case c_fx:
        case c_fy:
        case c_fz:
            return atom.properties->forces[col-5];
        }
    }else if(role == Qt::CheckStateRole){
        int col = colMap[index.column()];
        const auto& atom = curStep->at(index.row());
        // return bool*2 so true becomes CheckState::Checked
        switch(col){
        case c_hidden:
            return 2*atom.properties->flags[AtomProperties::Hidden];
        case c_fixx:
            return 2*atom.properties->flags[AtomProperties::FixX];
        case c_fixy:
            return 2*atom.properties->flags[AtomProperties::FixY];
        case c_fixz:
            return 2*atom.properties->flags[AtomProperties::FixZ];
        }
    }
    return QVariant{};
}

bool AtomModel::setData(const QModelIndex &index, const QVariant &value, int role)
{
    if (data(index, role) != value) {
        if(role == Qt::EditRole){
            int col = colMap[index.column()];
            auto atom = curStep->at(index.row());
            switch(col){
            case c_name:
                vApp.invokeOnStep([](Step &s, int i, const std::string &name){
                    s.at(i).name = name;
                }, index.row(), value.toString().toStdString());
                break;
            case c_x:
            case c_y:
            case c_z:
                // TODO: offload to vApp
                atom.coord[col-1] = value.toDouble();
                emit vApp.stepChanged(*vApp.curStep);
                break;
            case c_charge:
                vApp.invokeOnStep([](Step &s, int i, double charge){
                    s.at(i).properties->charge = charge;
                }, index.row(), value.toDouble());
                break;
            case c_fx:
            case c_fy:
            case c_fz:
                vApp.invokeOnStep([](Step &s, int i, int dir, double charge){
                    s.at(i).properties->forces[dir] = charge;
                }, index.row(), col-c_fx, value.toDouble());
                break;
            }
        }else if(role == Qt::CheckStateRole){
            int col = colMap[index.column()];
            auto atom = curStep->at(index.row());
            switch(col){
            case c_hidden:
                vApp.invokeOnStep([](Step &s, int i, bool val){
                    s.at(i).properties->flags[AtomProperties::Hidden] = val;
                }, index.row(), value == Qt::Checked);
                break;
            case c_fixx:
                vApp.invokeOnStep([](Step &s, int i, bool val){
                    s.at(i).properties->flags[AtomProperties::FixX] = val;
                }, index.row(), value == Qt::Checked);
                break;
            case c_fixy:
                vApp.invokeOnStep([](Step &s, int i, bool val){
                    s.at(i).properties->flags[AtomProperties::FixY] = val;
                }, index.row(), value == Qt::Checked);
                break;
            case c_fixz:
                vApp.invokeOnStep([](Step &s, int i, bool val){
                    s.at(i).properties->flags[AtomProperties::FixZ] = val;
                }, index.row(), value == Qt::Checked);
                break;
            }
        }else{
            return false;
        }
        emit dataChanged(index, index, QVector<int>() << role);
        return true;
    }
    return false;
}

Qt::ItemFlags AtomModel::flags(const QModelIndex &index) const
{
    if (!index.isValid() || !curStep)
        return Qt::NoItemFlags;

    if(colMap[index.column()] >= 8)
        return Qt::ItemIsEnabled | Qt::ItemIsSelectable | Qt::ItemIsUserCheckable;
    return Qt::ItemIsEditable | Qt::ItemIsEnabled | Qt::ItemIsSelectable;
}
