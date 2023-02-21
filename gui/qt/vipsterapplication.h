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

// Singleton to manage process state
class Application: public QObject
{
    Q_OBJECT

public:
    /* Accessor function
     * returns single program-wide instance
     * static instance is created before first access, has life-time of whole program
     */
    static Application& instance() {
        static Application app{};
        return app;
    }

    // TODO: use weak_ptr for dependent state data

    // expose config read from file
private:
public:
    ConfigState config{};
signals:
    void configChanged(Vipster::ConfigState &cfg);

    // Molecule state
private:
public:
    std::list<Molecule> molecules{};
    MolModel molModel{};
    Molecule *curMol{nullptr};
public:
    void newMol(Molecule &&mol);
signals:
    void molListChanged(const std::list<Vipster::Molecule> &molecules);
    void molChanged(int idx, Vipster::Molecule &mol);
    void activeMolChanged(Vipster::Molecule &mol);

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
public:
    template<class S, class F, class... Args>
    std::invoke_result_t<F&&, S&, Args&&...> invokeImpl(S &s, F &&f, Args &&...args)
    {
        if constexpr (!std::is_void_v<std::invoke_result_t<F&&, S&, Args&&...>>) {
            auto tmp = std::invoke(f, s, args...);
            emit stepChanged(*curStep);
            return tmp;
        } else {
            std::invoke(f, s, args...);
            emit stepChanged(*curStep);
        }
    }
public:
    Step *curStep{nullptr};
    Step::selection *curSel{nullptr};
    std::unique_ptr<Step::selection> copyBuf{}; // TODO: are lifetime semantics of selection sufficient??
    void setActiveStep(Step &step, Step::selection &sel);
    void selectionToCopy();

    template<class F, class... Args>
    std::invoke_result_t<F&&, Step&, Args&&...> invokeOnStep(F &&f, Args &&...args)
    {
        return invokeImpl(*curStep, std::forward<F>(f), std::forward<Args>(args)...);
    }
    template<class F, class... Args>
    std::invoke_result_t<F&&, Step::selection&, Args&&...> invokeOnSel(F &&f, Args &&...args)
    {
        return invokeImpl(*curSel, std::forward<F>(f), std::forward<Args>(args)...);
    }
signals:
    void activeStepChanged(Vipster::Step &step, Vipster::Step::selection &sel);
    void stepChanged(Step &step);
    void selChanged(Step::selection &sel);
    void copyBufChanged(Step::selection &buf);

public:
    // Data related to a loaded Step
    struct StepState{
        bool automatic_bonds{true};
        std::map<std::string,
                 std::tuple<Step::selection,
                            SelectionFilter,
                            std::shared_ptr<GUI::SelData>>> definitions;
    };
    std::map<Step*, StepState> stepdata{};

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
