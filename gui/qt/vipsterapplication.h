#ifndef VIPSTERAPPLICATION_H
#define VIPSTERAPPLICATION_H

#include <list>
#include <map>
#include <type_traits>
#include <QAction>
#include <QObject>

#include "seldata.h"
#include "molmodel.h"

#include "vipster/configfile.h"
#include "vipster/molecule.h"
#include "vipster/parameters.h"
#include "vipster/presets.h"

#define vApp Vipster::Application::instance()

namespace Vipster{

/* Roadmap:
 * - encapsulates application-wide data and state
 *   - data example: loaded molecules
 *   - state example: currently selected step
 * - widgets only retain state that is required and specific for their own operation
 *   - references to application state inside of widgets should be as temporary as possible
 * - state is opaquely encapsulated here
 *   - data-members should be private
 *   - setters and getters are provided as required
 * - Qt-signals are preferrable over getters
 *   -> setter shall emit a signal, interested widgets connect to that signal
 * - signals should be as atomic as possible, e.g. distinguish between:
 *   - modification of a step
 *   - selecting a different step
 *   - providing a new step
 */
// TODO: hide state (curMol, curStep etc.) behind invoke interface
// TODO: provide const ref only via signals

// Singleton to manage process state
class Application: public QObject
{
    Q_OBJECT
private:
    // Functional abstraction for managed state
    // TODO: provide undo-mechanism via this function
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

public:
    /* Accessor function
     * returns single program-wide instance
     * static instance is created before first access, has life-time of whole program
     */
    static Application& instance() {
        static Application app{};
        return app;
    }

    // expose config read from file
    // TODO: split into multiple parts -> better signalling
private:
    ConfigState config_val{};
public:
    const ConfigState &config();
    template<typename F, typename... Args>
    auto invokeOnConfig(F &&f, Args &&...args)
    {
        if constexpr (!std::is_void_v<decltype(invokeImpl(config_val, std::forward<F>(f), std::forward<Args>(args)...))>) {
            auto tmp = invokeImpl(config_val, std::forward<F>(f), std::forward<Args>(args)...);
            emit configChanged(config_val);
            return tmp;
        } else {
            invokeImpl(config_val, std::forward<F>(f), std::forward<Args>(args)...);
            emit configChanged(config_val);
        }
    }
signals:
    void configChanged(const Vipster::ConfigState &cfg);

    // Molecule state
private:
public:
    std::list<Molecule> molecules{};
    MolModel molModel{}; // TODO: what should this be used for?
    Molecule &curMol(); // TODO: return const
    void setActiveMol(Molecule &m);
private:
    Molecule *pCurMol{nullptr};
public:
    void newMol(Molecule &&mol);
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
    template<typename F, typename... Args>
    void invokeOnTrajec(F &&f, Args &&...args)
    {
        static_assert(std::is_void_v<decltype(invokeImpl(*pCurStep, std::forward<F>(f), std::forward<Args>(args)...))>);
        for (auto &step: pCurMol->getSteps()) {
            invokeImpl(step, std::forward<F>(f), std::forward<Args>(args)...);
            emit stepChanged(step);
        }
    }
signals:
    void molListChanged(const std::list<Vipster::Molecule> &molecules);
    void molChanged(const Vipster::Molecule &mol);
    void activeMolChanged(const Vipster::Molecule &mol);

    // Additional data
private:
public:
    std::list<Preset> presets{};
    std::list<Parameter> parameters{};
    std::list<std::unique_ptr<const BaseData>> data{};
public:
    void newPreset(Preset &&p);
    void newParameter(Parameter &&p);
    void newData(std::unique_ptr<const Vipster::BaseData> &&d);
    void newIOData(IOTuple &&t);
signals:
    void presetListChanged(void);
    void parameterListChanged(void);
    void dataListChanged(void);

    // Currently active state
private:
    Step *pCurStep{nullptr};
    Step::selection *pCurSel{nullptr};
public:
    const Step& curStep();
    const Step::selection& curSel();
    void setActiveStep(Step &step, Step::selection &sel); // TODO: prefer an index based interface?

    // TODO: integrate with OS-specific copy-paste functionality instead
    void selectionToCopy();
    Step copyBuf{};

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
    template<typename F, typename... Args>
    auto invokeOnSel(F &&f, Args &&...args)
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
    void updateSelection(SelectionFilter filter);
signals:
    void activeStepChanged(const Vipster::Step &step, const Vipster::Step::selection &sel);
    void stepChanged(const Vipster::Step &step);
    void selChanged(const Vipster::Step::selection &sel);
    void copyBufChanged(const Vipster::Step &buf);

public:
    // GUI-Data related to a loaded Step
    struct StepState{
        bool automatic_bonds;
        Step::formatter formatter;
        std::map<QWidget *, std::tuple<Step::selection, std::shared_ptr<GUI::SelData>>> selections;
        std::map<std::string,
                 std::tuple<Step::const_selection,
                            SelectionFilter,
                            std::shared_ptr<GUI::SelData>>> definitions;
    };
    // TODO: should return const ref
    StepState& getState(const Step &s);
private:
    // TODO: use weak_ptr
    std::map<const Step*, StepState> stepdata{};

private:
    // don't allow user-side construction
    Application();
    Application(const Application&) = delete;
    Application(Application&&) = delete;
    Application& operator=(const Application&) = delete;
    Application& operator=(Application&&) = delete;
};

}

#endif // VIPSTERAPPLICATION_H
