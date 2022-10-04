#ifndef LAMMPSWIDGET_LMP_H
#define LAMMPSWIDGET_LMP_H

#include "../basewidget.h"
#include "lammpswidget_aux/forcefield.lmp.h"

#include <filesystem>

namespace Ui {
class LammpsWidget;
}

class LammpsWidget : public BaseWidget
{
    Q_OBJECT

public:
    explicit LammpsWidget(QWidget *parent = nullptr);
    ~LammpsWidget();

private slots:
    void on_runButton_clicked();
    void on_helpButton_clicked();

    void on_ffPrepare_clicked();

    void on_ffSel_currentIndexChanged(int index);

private:
    Ui::LammpsWidget *ui;
    ForceFields forcefields;
    void mkScript(const Vipster::Step &curStep, const ForceField &FF, const std::filesystem::path &tempdir);
    void mkGeom(const Vipster::Step &curStep, const ForceField &FF, const std::filesystem::path &tempdir);
};

#endif // LAMMPSWIDGET_LMP_H
