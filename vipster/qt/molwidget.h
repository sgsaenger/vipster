#ifndef MOLWIDGET_H
#define MOLWIDGET_H

#include <QWidget>
#include "basewidget.h"
#include "molecule.h"

namespace Ui {
class MolWidget;
}

class MolWidget : public BaseWidget
{
    Q_OBJECT

public:
    explicit MolWidget(QWidget *parent = nullptr);
    ~MolWidget() override;
    void updateWidget(Vipster::guiChange_t change) override;
    void registerMol(const std::string& name);
    Vipster::AtomFmt getAtomFmt();
    Vipster::CdmFmt getCellFmt();

private slots:
    void on_cellEnabledBox_toggled(bool checked);
    void on_cellFmt_currentIndexChanged(int idx);
    void on_cellDimBox_valueChanged(double cdm);
    void on_cellVecTable_cellChanged(int row, int column);
    void on_cellTrajecButton_clicked();
    void on_atomTable_cellChanged(int row, int column);
    void on_atomFmtBox_currentIndexChanged(int index);
    void on_atomFmtButton_clicked();
    void on_atomTableButton_toggled(bool checked);
    void on_molList_currentIndexChanged(int index);
    void on_atomTable_itemSelectionChanged();
    void on_kFmtButton_clicked();
    void on_bands_stateChanged(int);
    void on_crystal_stateChanged(int);

    void mpg_change();

    void on_discretetable_itemSelectionChanged();
    void on_actionNew_K_Point_triggered();
    void on_actionDelete_K_Point_triggered();

    void on_discretetable_cellChanged(int row, int column);

private:
    void fillAtomTable(void);
    void fillCell(void);
    void fillKPoints(void);
    void setSelection(void);
    Ui::MolWidget *ui;
    Vipster::Step curStep;
    Vipster::Molecule* curMol;
    int curKPoint{-1};
    bool atomsOutdated{true};
    static constexpr const char* inactiveKpoints[] = {"Gamma", "Monkhorst-Pack grid", "Discrete"};
    static constexpr const char* activeKpoints[] = {"Gamma (active)",
                                                    "Monkhorst-Pack grid (active)",
                                                    "Discrete (active)"};
    static constexpr const char* inactiveFmt[] = {"Bohr", "Angstrom", "Crystal", "Alat"};
    static constexpr const char* activeFmt[] = {"Bohr (active)",
                                                "Angstrom (active)",
                                                "Crystal (active)",
                                                "Alat (active)"};
};

#endif // MOLWIDGET_H
