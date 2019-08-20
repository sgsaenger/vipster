#ifndef BONDMODEL_H
#define BONDMODEL_H

#include <QAbstractTableModel>
#include "step.h"

class BondWidget;
class BondModel: public QAbstractTableModel
{
    Q_OBJECT
public:
    explicit BondModel(BondWidget* parent=nullptr);
    void setStep(Vipster::Step* curStep);

    // Header:
    QVariant headerData(int section, Qt::Orientation orientation,
                        int role = Qt::DisplayRole) const override;

    // Basic functionality:
    int rowCount(const QModelIndex &parent = QModelIndex()) const override;
    int columnCount(const QModelIndex &parent = QModelIndex()) const override;

    QVariant data(const QModelIndex &index, int role = Qt::DisplayRole) const override;

    // Editable:
    bool setData(const QModelIndex &index, const QVariant &value,
                 int role = Qt::EditRole) override;

    Qt::ItemFlags flags(const QModelIndex& index) const override;
private:
    BondWidget *parent;
    Vipster::Step *curStep{nullptr};
    const std::vector<Vipster::Bond> *curBonds{nullptr};
    QStringList colNames = {"Atoms", "Length", "Type", "Color"};
};

#endif // BONDMODEL_H
