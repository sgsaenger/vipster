#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QAbstractButton>
#include <QDir>
#include <QTimer>
#include <vector>
#include "io.h"
#include "stepsel.h"
#include "configfile.h"
#include "../common/guiwrapper.h"
#include "mainwidgets/paramwidget.h"
#include "mainwidgets/configwidget.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QString path, Vipster::ConfigState& state,
                        std::vector<Vipster::IO::Data> &&d={},
                        QWidget *parent = nullptr);
    ~MainWindow() override;
    // Molecule and Step data
    Vipster::Molecule* curMol{nullptr};
    Vipster::Step* curStep{nullptr};
    Vipster::Step::selection* curSel{nullptr};
    Vipster::Step copyBuf{};
    void updateWidgets(Vipster::GUI::change_t change);
    void newData(Vipster::IO::Data&& d);
    void setMultEnabled(bool);
    struct MolExtras{
        int curStep{-1};
        Vipster::GUI::PBCVec mult{1,1,1};
    };
    struct StepExtras{
        std::unique_ptr<Vipster::Step::selection> sel{nullptr};
        std::map<std::string, Vipster::Step::selection> def{};
    };
    std::list<Vipster::Molecule> molecules;
    std::map<Vipster::Molecule*, MolExtras> moldata;
    std::map<Vipster::Step*, StepExtras> stepdata;
    // Parameter data
    std::map<Vipster::IOFmt, QMenu*> paramMenus;
    ParamWidget* paramWidget;
    const decltype (ParamWidget::params)& getParams() const noexcept;
    // Config data
    std::map<Vipster::IOFmt, QMenu*> configMenus;
    ConfigWidget* configWidget;
    const decltype (ConfigWidget::configs)& getConfigs() const noexcept;
    // Extra data
    std::list<std::unique_ptr<const Vipster::BaseData>> data;
    // expose configstate read from file
    Vipster::ConfigState    &state;
    Vipster::PeriodicTable  &pte;
    Vipster::Settings       &settings;
    Vipster::IO::Parameters &params;
    Vipster::IO::Configs    &configs;
    // GL helpers for additional render-data
    void addExtraData(Vipster::GUI::Data* dat);
    void delExtraData(Vipster::GUI::Data* dat);
    const Vipster::GUI::GlobalData& getGLGlobals();

public slots:
    void setMol(int i);
    void setStep(int i);
    void setMult(int i);
    void setBondMode(int i);
    void stepBut(QAbstractButton *but);
    void about();
    void newMol();
    void newMol(QAction *sender);
    void loadMol();
    void saveMol();
    void editAtoms(QAction *sender);
    void loadParam();
    void saveParam();
    void loadConfig();
    void saveConfig();
    void saveScreenshot();

private slots:
    void on_molList_currentIndexChanged(int index);

private:
    void setupUI(void);
    void registerMol(const std::string& name);

    Ui::MainWindow *ui;
    QDir path{};
    QTimer playTimer{};
    std::vector<BaseWidget*> baseWidgets;
    std::vector<BaseWidget*> toolWidgets;
};
#endif // MAINWINDOW_H
