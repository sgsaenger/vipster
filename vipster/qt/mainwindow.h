#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QAbstractButton>
#include <vector>
#include "iowrapper.h"
#include "../common/guiwrapper.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
    explicit MainWindow(Vipster::IO::Data &&d, QWidget *parent = nullptr);
    ~MainWindow();
    Vipster::Molecule *curMol{nullptr};
    Vipster::StepProper *curStep{nullptr};
    Vipster::IO::BaseParam *curParam{nullptr};
    Vipster::AtomFmt getFmt();
    void setFmt(int i, bool apply, bool scale);
    void updateWidgets(Vipster::Change change);
    void newData(Vipster::IO::Data&& d);

public slots:
    void setMol(int i);
    void setStep(int i);
    void setParam(int i);
    void stepBut(QAbstractButton *but);
    void about(void);
    void newMol();
    void loadMol();
    void editAtoms(void);

private:
    void setupUI(void);
    void newParam(Vipster::IOFmt fmt, std::unique_ptr<Vipster::IO::BaseParam> &&p);
    void newMol(Vipster::Molecule &&m);

    Vipster::AtomFmt fmt{Vipster::AtomFmt::Bohr};
    Ui::MainWindow *ui;
    std::vector<Vipster::Molecule> molecules;
    std::vector<std::pair<Vipster::IOFmt,
            std::unique_ptr<Vipster::IO::BaseParam>>> params;
    std::vector<std::pair<Vipster::IOFmt,
            std::unique_ptr<Vipster::IO::BaseConfig>>> configs;
};

#endif // MAINWINDOW_H
