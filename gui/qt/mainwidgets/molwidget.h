#ifndef MOLWIDGET_H
#define MOLWIDGET_H

#include <QWidget>
#include <QItemSelection>
#include "vipster/molecule.h"
#include "../basewidget.h"
#include "molwidget_aux/atommodel.h"
#include "molwidget_aux/bondmodel.h"
#include "molwidget_aux/cellmodel.h"

namespace Ui {
class MolWidget;
}

class MolWidget : public BaseWidget
{
    Q_OBJECT

public:
    explicit MolWidget(QWidget *parent = nullptr);
    ~MolWidget() override;
    void updateWidget(Vipster::GUI::change_t change) override;

private slots:
    // atom slots
    void on_atomFmtBox_currentIndexChanged(int index);
    void on_atomFmtButton_clicked();
    void atomSelectionChanged(const QItemSelection &selected, const QItemSelection &deselected);
    void on_atomHelpButton_clicked();

    // cell slots
    void on_cellEnabledBox_toggled(bool checked);
    void on_cellFmt_currentIndexChanged(int idx);
    void on_cellDimBox_valueChanged(double cdm);
    void on_cellTrajecButton_clicked();

    // kpoint slots
    void on_kFmtButton_clicked();
    void on_bands_stateChanged(int);
    void on_crystal_stateChanged(int);
    void mpg_change();
    void on_actionNew_K_Point_triggered();
    void on_actionDelete_K_Point_triggered();
    void on_discretetable_itemSelectionChanged();
    void on_discretetable_cellChanged(int row, int column);

    // bond slots
    void on_bondSetButton_clicked();
    void on_bondHelpButton_clicked();
    void on_bondModeBox_currentIndexChanged(int index);
    void on_ovlpTable_itemSelectionChanged();

    // pte slots
    void on_clearTableButton_clicked();
    void on_newElemButton_clicked();

private:
    void setActiveStep(Vipster::Step &step, Vipster::Step::selection &sel);
    void updateStep(Vipster::Step &step);
    void updateSelection(Vipster::Step::selection &sel);

    void setActiveMol(Vipster::Molecule &mol);

    bool scale();
    void checkOverlap(void);
    void fillCell(void);
    void fillKPoints(void);
    Ui::MolWidget *ui;
    Vipster::Step *curStep;
    Vipster::Step::selection *curSel;
    std::unique_ptr<Vipster::Step::formatter> ownStep;
    Vipster::Molecule* curMol;
    AtomModel atomModel{this};
    BondModel bondModel{this};
    friend class CellModel;
    CellModel cellModel{this};
    QList<QAction*> headerActions;
    int curKPoint{-1};
    static constexpr const char* inactiveKpoints[] = {"Gamma", "Monkhorst-Pack grid", "Discrete"};
    static constexpr const char* activeKpoints[] = {"Gamma (active)",
                                                    "Monkhorst-Pack grid (active)",
                                                    "Discrete (active)"};
    static constexpr const char* inactiveFmt[] = { "Crystal", "Alat", "Angstrom", "Bohr" };
    static constexpr const char* activeFmt[] = {"Crystal (active)",
                                                "Alat (active)",
                                                "Angstrom (active)",
                                                "Bohr (active)"};
};

#endif // MOLWIDGET_H
