#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QDir>
#include <QSplitter>
#include <vector>

#include "fileio.h"
#include "step.h"
#include "configfile.h"

#include "viewport.h"
#include "guiwrapper.h"
#include "guidata.h"
#include "mainwidgets/paramwidget.h"
#include "mainwidgets/presetwidget.h"

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
    // TODO: store StepData with weak_ptr
    struct StepState{
        std::map<std::string,
                 std::tuple<Vipster::Step::selection,
                            Vipster::SelectionFilter,
                            std::shared_ptr<Vipster::GUI::SelData>>> definitions;
        bool automatic_bonds{true};
    };
    std::map<Vipster::Step*, StepState> stepdata{};
    Vipster::Molecule* curMol{nullptr};
    Vipster::Step* curStep{nullptr};
    Vipster::Step::selection* curSel{nullptr};
    std::unique_ptr<Vipster::Step::selection> copyBuf{};
    void updateWidgets(Vipster::GUI::change_t change);
    void newData(Vipster::IO::Data&& d);
    // Parameter data
    std::map<const Vipster::IO::Plugin*, QMenu*> paramMenus;
    ParamWidget* paramWidget;
    const decltype (ParamWidget::params)& getParams() const noexcept;
    // Preset data
    std::map<const Vipster::IO::Plugin*, QMenu*> presetMenus;
    PresetWidget* presetWidget;
    const decltype (PresetWidget::presets)& getPresets() const noexcept;
    // Extra data
    std::list<std::unique_ptr<const Vipster::BaseData>> data;
    // expose configstate read from file
    Vipster::ConfigState    &state;
    Vipster::PeriodicTable  &pte;
    Vipster::Settings       &settings;
    Vipster::IO::Plugins    &plugins;
    Vipster::IO::Parameters &params;
    Vipster::IO::Presets    &presets;

public slots:
    void about();
    void newMol();
    void newMol(QAction *sender);
    void loadMol();
    void saveMol();
    void editAtoms(QAction *sender);
    void loadParam();
    void saveParam();
    void loadPreset();
    void savePreset();
    void saveScreenshot();

private:
    void setupUI(void);
    void registerMol(const std::string& name);

    Ui::MainWindow *ui;
    QDir path{};
    QSplitter *vsplit;
    std::vector<QSplitter*> hsplits;
    std::vector<BaseWidget*> mainWidgets;
    std::vector<BaseWidget*> toolWidgets;
};
#endif // MAINWINDOW_H
