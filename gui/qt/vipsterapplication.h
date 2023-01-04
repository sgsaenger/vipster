#ifndef VIPSTERAPPLICATION_H
#define VIPSTERAPPLICATION_H

#include <list>
#include <map>
#include "seldata.h"
#include "vipster/configfile.h"
#include "vipster/molecule.h"
#include "vipster/parameters.h"
#include "vipster/presets.h"

#define vApp Vipster::GUI::Application::instance()

namespace Vipster::GUI{

// Singleton to manage process state
class Application
{
public:
    static Application& instance() {
        static Application app{};
        return app;
    }

    // TODO: use weak_ptr for dependent state data

    // expose config read from file
    Vipster::ConfigState config{};

    // Currently open data
    std::list<Vipster::Molecule> molecules{};
    std::list<Vipster::Preset> presets{};
    std::list<Vipster::Parameter> parameters{};
    std::list<std::unique_ptr<const Vipster::BaseData>> data{};
    std::unique_ptr<Vipster::Step::selection> copyBuf{};

    // Currently active state
    Vipster::Molecule *curMol{nullptr};
    Vipster::Step *curStep{nullptr};
    Vipster::Step::selection *curSel{nullptr};

    // Data related to a loaded Step
    struct StepState{
        bool automatic_bonds{true};
        std::map<std::string,
                 std::tuple<Step::selection,
                            SelectionFilter,
                            std::shared_ptr<SelData>>> definitions;
    };
    std::map<Vipster::Step*, StepState> stepdata{};

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
