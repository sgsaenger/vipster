#ifndef CELLMODEL_H
#define CELLMODEL_H

#include <QAbstractTableModel>
#include "vipster/step.h"

class MolWidget;
class CellModel : public QAbstractTableModel
{
    Q_OBJECT
public:
    explicit CellModel(MolWidget *parent = nullptr);
    void setStep(Vipster::Step *curStep);

    Vipster::Mat getVec();

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

    Qt::ItemFlags flags(const QModelIndex &index) const override;
private:
    MolWidget *parent;
    Vipster::Step *curStep{nullptr};
    Vipster::Mat vec{};
};

#endif // CELLMODEL_H
