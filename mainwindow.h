#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QMessageBox>
#include <vector>
#include <array>
#include "libvipster.h"

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
    inline Vipster::Molecule& curMol();

public slots:
    void setMol(int i);
    void setStep(int i=-1);
    void about(void);
    void newMol(Vipster::Molecule m=Vipster::Molecule());
    void editAtoms(void);

private slots:
    void on_cellDimBox_valueChanged(double arg1);
    void on_cellVecTable_cellChanged(int row, int column);
    void on_atomTable_cellChanged(int row, int column);

private:
    Ui::MainWindow *ui;
    std::vector<Vipster::Molecule> molecules;
    uint molIdx;
};

#endif // MAINWINDOW_H
