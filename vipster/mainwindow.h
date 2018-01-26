#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QAbstractButton>
#include <vector>
#include "molecule.h"

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

public slots:
    void setMol(void);
    void setMol(int i);
    void setStep(void);
    void setStep(int i);
    void stepBut(QAbstractButton *but);
    void about(void);
    void newMol();
    void newMol(const Vipster::Molecule &m);
    void newMol(Vipster::Molecule &&m);
    void editAtoms(void);

private:
    void setupUI(void);

    Ui::MainWindow *ui;
    std::vector<Vipster::Molecule> molecules;
    Vipster::Molecule *curMol = nullptr;
    Vipster::StepProper *curStep = nullptr;
};

#endif // MAINWINDOW_H
