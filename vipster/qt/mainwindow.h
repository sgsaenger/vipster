#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QAbstractButton>
#include <QDir>
#include <vector>
#include "iowrapper.h"
#include "stepsel.h"
#include "../common/guiwrapper.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
    explicit MainWindow(std::vector<Vipster::IO::Data> &&d, QWidget *parent = nullptr);
    ~MainWindow() override;
    Vipster::Molecule* curMol{nullptr};
    Vipster::StepProper* curStep{nullptr};
    Vipster::BaseParam* curParam{nullptr};
    Vipster::StepSelection* curSel{nullptr};
    Vipster::AtomFmt getFmt();
    void setFmt(int i, bool apply, bool scale);
    void updateWidgets(uint8_t change);
    void newData(Vipster::IO::Data&& d);
    std::vector<std::pair< Vipster::IOFmt, std::unique_ptr<Vipster::BaseParam>>> params;
    std::vector<std::pair< Vipster::IOFmt, std::unique_ptr<Vipster::BaseConfig>>> configs;
    std::vector<Vipster::Molecule> molecules;

public slots:
    void setMol(int i);
    void setStep(int i);
    void setParam(int i);
    void stepBut(QAbstractButton *but);
    void about(void);
    void newMol();
    void loadMol();
    void saveMol();
    void editAtoms(void);
    void loadParam();

private:
    void setupUI(void);

    Vipster::AtomFmt fmt{Vipster::AtomFmt::Bohr};
    Ui::MainWindow *ui;
    QDir path{QString{std::getenv("HOME")}};
};

class BaseWidget{
public:
    BaseWidget();
    void triggerUpdate(uint8_t change);
    virtual void updateWidget(uint8_t){}
    virtual ~BaseWidget()=default;
protected:
    bool updateTriggered{false};
    MainWindow* master;
};

#endif // MAINWINDOW_H
