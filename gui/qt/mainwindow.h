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
#include "molmodel.h"


/* Main window and application state for Qt-based GUI
 *
 * Singleton class - access globally via `vApp` macro
 * - encapsulates application-wide data and state
 * - synchronizes state across GUI via signals & slots
 * - sub-widgets only retain state that is required for their own operation
 *   -> connect to relevant signals
 */
#define vApp MainWindow::instance()
class MainWindow : public QMainWindow
{
    Q_OBJECT

    // Singleton: exactly one instance can/must exist
public:
    explicit MainWindow(Vipster::ConfigState state, QDir path);
    static MainWindow &instance() {
        return *MainWindow::_instance;
    }
private:
    static MainWindow *_instance;
    // Singleton: no copy/move
    MainWindow(const MainWindow&) = delete;
    MainWindow& operator=(const MainWindow&) = delete;
    MainWindow(MainWindow&&) = delete;
    MainWindow& operator=(MainWindow&&) = delete;


    /* Global application configuration
     */
private:
    Vipster::ConfigState _config{};
public:
    const Vipster::ConfigState &config();
    template<typename F, typename... Args>
    auto invokeOnConfig(F &&f, Args &&...args)
    {
        if constexpr (!std::is_void_v<decltype(invokeImpl(_config, std::forward<F>(f), std::forward<Args>(args)...))>) {
            auto tmp = invokeImpl(_config, std::forward<F>(f), std::forward<Args>(args)...);
            emit configChanged(_config);
            return tmp;
        } else {
            invokeImpl(_config, std::forward<F>(f), std::forward<Args>(args)...);
            emit configChanged(_config);
        }
    }
signals:
    void configChanged(const Vipster::ConfigState &cfg);


    /* File operations
     */
public:
    void newIOData(Vipster::IOTuple &&t);
private:
    void readFile();
    void saveFile();


    /* Molecule & Step data
     */
private:
public:
    std::list<Vipster::Molecule> molecules{};
    MolModel molModel{};                        // TODO: what should this be used for?
    Vipster::Molecule &curMol();                // TODO: return const
    const Vipster::Step& curStep();
    const Vipster::Step::selection& curSel();
    // Update currently active state, to be called by active viewport only
    void setActiveMol(Vipster::Molecule &m);
    void setActiveStep(Vipster::Step &step, Vipster::Step::selection &sel); // TODO: prefer an index based interface?
    // TODO: integrate with OS-specific copy-paste functionality instead
    void selectionToCopy();
    Vipster::Step copyBuf{};
private:
    Vipster::Molecule        * pCurMol{};
    Vipster::Step            * pCurStep{};
    Vipster::Step::selection * pCurSel{};
public:
    // Load new molecule
    void newMol(Vipster::Molecule &&mol);
    // Modify molecule properties (e.g. PTE, KPoints)
    template<typename F, typename... Args>
    auto invokeOnMol(F &&f, Args &&...args)
    {
        if constexpr (!std::is_void_v<decltype(invokeImpl(*pCurMol, std::forward<F>(f), std::forward<Args>(args)...))>) {
            auto tmp = invokeImpl(*pCurMol, std::forward<F>(f), std::forward<Args>(args)...);
            emit molChanged(*pCurMol);
            return tmp;
        } else {
            invokeImpl(*pCurMol, std::forward<F>(f), std::forward<Args>(args)...);
            emit molChanged(*pCurMol);
        }
    }
    // Modify currently active Step
    template<typename F, typename... Args>
    auto invokeOnStep(F &&f, Args &&...args)
    {
        if constexpr (!std::is_void_v<decltype(invokeImpl(*pCurStep, std::forward<F>(f), std::forward<Args>(args)...))>) {
            auto tmp = invokeImpl(*pCurStep, std::forward<F>(f), std::forward<Args>(args)...);
            emit stepChanged(*pCurStep);
            return tmp;
        } else {
            invokeImpl(*pCurStep, std::forward<F>(f), std::forward<Args>(args)...);
            emit stepChanged(*pCurStep);
        }
    }
    // Modify selected subset of active Step
    template<typename F, typename... Args>
    auto invokeOnSelection(F &&f, Args &&...args)
    {
        if constexpr (!std::is_void_v<decltype(invokeImpl(*pCurSel, std::forward<F>(f), std::forward<Args>(args)...))>) {
            auto tmp = invokeImpl(*pCurSel, std::forward<F>(f), std::forward<Args>(args)...);
            emit stepChanged(*pCurStep);
            return tmp;
        } else {
            invokeImpl(*pCurSel, std::forward<F>(f), std::forward<Args>(args)...);
            emit stepChanged(*pCurStep);
        }
    }
    // Modify all Steps of active Mol
    template<typename F, typename... Args>
    void invokeOnTrajec(F &&f, Args &&...args)
    {
        static_assert(std::is_void_v<decltype(invokeImpl(*pCurStep, std::forward<F>(f), std::forward<Args>(args)...))>);
        for (auto &step: pCurMol->getSteps()) {
            invokeImpl(step, std::forward<F>(f), std::forward<Args>(args)...);
            emit stepChanged(step);
        }
    }
    // select a new subset of active Step
    void updateSelection(Vipster::SelectionFilter filter);
signals:
    void molListChanged(const std::list<Vipster::Molecule> &molecules);
    void molChanged(const Vipster::Molecule &mol);
    void activeMolChanged(const Vipster::Molecule &mol);
    void activeStepChanged(const Vipster::Step &step, const Vipster::Step::selection &sel);
    void stepChanged(const Vipster::Step &step);
    void selChanged(const Vipster::Step::selection &sel);
    void copyBufChanged(const Vipster::Step &buf);
public:
    // GUI-Data related to a loaded Step
    struct StepState{
        bool automatic_bonds;
        Vipster::Step::formatter formatter;
        std::map<ViewPort*, Vipster::Step::selection> selections;
        std::map<std::string,
                 std::tuple<Vipster::Step::const_selection,
                            Vipster::SelectionFilter,
                            std::shared_ptr<Vipster::GUI::SelData>>> definitions;
    };
    // TODO: should return const ref
    StepState& getState(const Vipster::Step &s);
private:
    // TODO: use weak_ptr or other mechanism to allow clean-up
    std::map<const Vipster::Step*, StepState> stepdata{};


    // Parameter data
public:
    const decltype (ParamWidget::params)& getParams() const noexcept;
private:
    std::map<const Vipster::Plugin*, QMenu*> paramMenus;
    ParamWidget* paramWidget;
    void loadParam();
    void saveParam();
signals:
    void parameterListChanged(void);

    // IO-Preset data
public:
    const decltype (PresetWidget::presets)& getPresets() const noexcept;
private:
    std::map<const Vipster::Plugin*, QMenu*> presetMenus;
    PresetWidget* presetWidget;
    void loadPreset();
    void savePreset();
signals:
    void presetListChanged(void);

    // Additional data (e.g. density data)
public:
// private:
    std::list<std::unique_ptr<const Vipster::BaseData>> data{};
signals:
    void dataListChanged(void);

    // Viewports
public:
    ViewPort *curVP{};
    std::vector<ViewPort*> viewports;
    void setActiveViewport(ViewPort* sender);
    void splitViewportHoriz(ViewPort* sender);
    void splitViewportVert(ViewPort* sender);
    void closeViewport(ViewPort* sender);

    // Screenshot
public:
    void saveScreenshot(QString fn);
    void saveScreenshot();
    void saveScreenshots();

    // Implementation details
private:
    QDir path{};

    QSplitter *vsplit;
    std::vector<QSplitter*> hsplits;

    // GUI instantiation helpers
    void setupMainWidgets();
    void setupToolWidgets();
    void setupFileMenu();
    void setupEditMenu();
    void setupHelpMenu();
    void setupViewports();

    /* State manipulation
     *
     * Any operations modifying the global state shall be invoked via any of the functions implemented here.
     * This shall ensure triggering the required signals.
     * TODO: provide undo-mechanism via these functions
     */
private:
    template<typename S, typename F, typename... Args>
    std::invoke_result_t<F&&, S&, Args&&...> invokeImpl(S &s, F &&f, Args &&...args)
    {
        // TODO: if not required unless other functionality is introduced
        if constexpr (!std::is_void_v<std::invoke_result_t<F&&, S&, Args&&...>>) {
            return std::invoke(std::forward<F>(f), s, std::forward<Args>(args)...);
        } else {
            std::invoke(std::forward<F>(f), s, std::forward<Args>(args)...);
        }
    }
};
#endif // MAINWINDOW_H
