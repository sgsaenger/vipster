#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QAbstractButton>
#include <vector>
#include "ioplugin.h"
#include "../common/guiwrapper.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
    explicit MainWindow(const Vipster::IO::BaseData &d, QWidget *parent = nullptr);
    explicit MainWindow(Vipster::IO::BaseData &&d, QWidget *parent = nullptr);
    ~MainWindow();
    Vipster::Molecule *curMol{nullptr};
    Vipster::StepProper *curStep{nullptr};
    Vipster::IO::BaseParam *curParam{nullptr};
    Vipster::AtomFmt getFmt();
    void setFmt(int i, bool apply, bool scale);
    void updateWidgets(Vipster::Change change);

public slots:
    void setMol(int i);
    void setStep(int i);
    void setParam(int i);
    void stepBut(QAbstractButton *but);
    void about(void);
    void newMol();
    void newMol(const Vipster::Molecule &m);
    void newMol(Vipster::Molecule &&m);
    void newParam(const std::unique_ptr<Vipster::IO::BaseParam> &p);
    void newParam(std::unique_ptr<Vipster::IO::BaseParam> &&p);
    void editAtoms(void);

private:
    void setupUI(void);

    Vipster::AtomFmt fmt{Vipster::AtomFmt::Bohr};
    Ui::MainWindow *ui;
    std::vector<Vipster::Molecule> molecules;
    std::vector<std::unique_ptr<Vipster::IO::BaseParam>> params;
    std::vector<std::unique_ptr<Vipster::IO::BaseConfig>> configs;
};

#endif // MAINWINDOW_H
