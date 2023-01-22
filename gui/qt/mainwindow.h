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

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QString path,
                        Vipster::ConfigState& state, // TODO: remove this
                        QWidget *parent = nullptr);
    ~MainWindow() override;

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
    void saveScreenshots();

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
