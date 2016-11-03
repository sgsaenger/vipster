#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <vector>
#include <array>
#include <molecule.h>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    explicit MainWindow(Vipster::Molecule m, QWidget *parent = 0);
    ~MainWindow();
    Vipster::Molecule *curMol;
    Vipster::Step *curStep;

public slots:
    void setMol(void);
    void setMol(int i);
    void setStep(void);
    void setStep(int i);
    void about(void);
    void newMol(Vipster::Molecule m=Vipster::Molecule());
    void editAtoms(void);

private:
    Ui::MainWindow *ui;
    std::vector<Vipster::Molecule> molecules;
};

#endif // MAINWINDOW_H
