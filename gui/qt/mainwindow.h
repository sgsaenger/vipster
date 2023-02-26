#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QDir>
#include <QSplitter>
#include <vector>

#include "vipster/molecule.h"
#include "vipster/configfile.h"

#include "viewport.h"
#include "mainwidgets/paramwidget.h"
#include "mainwidgets/presetwidget.h"

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QString path);

    // Viewports
    std::vector<ViewPort*> viewports;
    enum VPChange{VP_CLOSE, VP_VSPLIT, VP_HSPLIT, VP_ACTIVE};
    void changeViewports(ViewPort* sender, VPChange change);
    ViewPort* curVP{nullptr};
    void updateWidgets(Vipster::GUI::change_t change);
    void newData(Vipster::IOTuple&& d);

    // Parameter data
    std::map<const Vipster::Plugin*, QMenu*> paramMenus;
    ParamWidget* paramWidget;
    const decltype (ParamWidget::params)& getParams() const noexcept;

    // Preset data
    std::map<const Vipster::Plugin*, QMenu*> presetMenus;
    PresetWidget* presetWidget;
    const decltype (PresetWidget::presets)& getPresets() const noexcept;

    void saveScreenshot(QString fn);

public slots:
    void loadMol();
    void saveMol();
    void loadParam();
    void saveParam();
    void loadPreset();
    void savePreset();
    void saveScreenshot();
    void saveScreenshots();

private:
    void setupMainWidgets();
    void setupToolWidgets();
    void setupFileMenu();
    void setupEditMenu();
    void setupHelpMenu();
    void setupViewports();

    QDir path{};
    QSplitter *vsplit;
    std::vector<QSplitter*> hsplits;
    std::vector<BaseWidget*> toolWidgets;
};
#endif // MAINWINDOW_H
