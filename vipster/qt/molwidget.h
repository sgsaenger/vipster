#ifndef MOLWIDGET_H
#define MOLWIDGET_H

#include <QWidget>
#include "mainwindow.h"

namespace Ui {
class MolWidget;
}

class MolWidget : public BaseWidget
{
    Q_OBJECT

public:
    explicit MolWidget(QWidget *parent = nullptr);
    ~MolWidget() override;
    void updateWidget(uint8_t change) override;
    void registerMol(const std::string& name);
    Vipster::AtomFmt getAtomFmt();
    Vipster::CdmFmt getCellFmt();

private slots:
    void on_cellEnabled_toggled(bool checked);
    void on_cellFmt_currentIndexChanged(int idx);
    void on_cellDimBox_valueChanged(double cdm);
    void on_cellVecTable_cellChanged(int row, int column);
    void on_atomTable_cellChanged(int row, int column);
    void on_atomFmtBox_currentIndexChanged(int index);
    void on_atomFmtButton_clicked();
    void on_atomTableButton_toggled(bool checked);
    void on_molList_currentIndexChanged(int index);
    void on_atomTable_itemSelectionChanged();

private:
    void fillAtomTable(void);
    void fillCell(void);
    void fillKPoints(void);
    void setSelection(void);
    Ui::MolWidget *ui;
    Vipster::Step *curStep;
    Vipster::Molecule* curMol;
    bool atomsOutdated{true};
};

#endif // MOLWIDGET_H
