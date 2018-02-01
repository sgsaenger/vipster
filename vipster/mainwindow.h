#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QAbstractButton>
#include <vector>
#include "molecule.h"
#include "guiwrapper.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    explicit MainWindow(const Vipster::Molecule &m, QWidget *parent = 0);
    explicit MainWindow(Vipster::Molecule &&m, QWidget *parent = 0);
    ~MainWindow();
    Vipster::Molecule *curMol{nullptr};
    Vipster::StepProper *curStep{nullptr};
    Vipster::AtomFmt getFmt();
    void setFmt(int i, bool apply, bool scale);
    void updateWidgets(Vipster::Change change);

private slots:
    void setMol(int i);
    void setStep(int i);
    void stepBut(QAbstractButton *but);
    void about(void);
    void newMol();
    void newMol(const Vipster::Molecule &m);
    void newMol(Vipster::Molecule &&m);
    void editAtoms(void);

private:
    void setupUI(void);

    Vipster::AtomFmt fmt{Vipster::AtomFmt::Bohr};
    Ui::MainWindow *ui;
    std::vector<Vipster::Molecule> molecules;
};

#endif // MAINWINDOW_H
