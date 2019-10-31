#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QDir>
#include <QSplitter>
#include <vector>

#include "io.h"
#include "stepsel.h"
#include "configfile.h"

#include "viewport.h"
#include "../common/guiwrapper.h"
#include "../common/guidata.h"
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
    // Viewports
    std::vector<ViewPort*> viewports;
    enum VPChange{VP_CLOSE, VP_VSPLIT, VP_HSPLIT, VP_ACTIVE};
    void changeViewports(ViewPort* sender, VPChange change);
    ViewPort* curVP{nullptr};
    Vipster::GUI::GlobalData globals{};
    // Molecule and Step data
    std::list<Vipster::Molecule> molecules;
    Vipster::Molecule* curMol{nullptr};
    Vipster::Step* curStep{nullptr};
    Vipster::Step::selection* curSel{nullptr};
    Vipster::Step copyBuf{};
    void updateWidgets(Vipster::GUI::change_t change);
    void newData(Vipster::IO::Data&& d);
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

public slots:
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

private:
    void setupUI(void);
    void registerMol(const std::string& name);

    Ui::MainWindow *ui;
    QDir path{};
    QSplitter *vsplit;
    std::vector<QSplitter*> hsplits;
    std::vector<BaseWidget*> baseWidgets;
    std::vector<BaseWidget*> toolWidgets;
};
#endif // MAINWINDOW_H
