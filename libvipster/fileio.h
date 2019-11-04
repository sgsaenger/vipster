#ifndef IOWRAPPER
#define IOWRAPPER

#include <map>
#include <optional>

#include "molecule.h"
#include "io/data.h"
#include "io/plugin.h"
#include "io/parameters.h"
#include "io/presets.h"

//TODO: check std::ios_base::sync_with_stdio(false)
namespace Vipster{
    // convenience containers and defaults
    namespace IO{
        using Plugins = std::vector<const Plugin*>;
        Plugins defaultPlugins();
        using Presets = std::map<const Plugin*, std::map<std::string, std::unique_ptr<BasePreset>>>;
        Presets defaultPresets(const Plugins& p);
        using Parameters = std::map<const Plugin*, std::map<std::string, std::unique_ptr<BaseParam>>>;
        Parameters defaultParams(const Plugins& p);
    }

    // read with format guess
    IO::Data readFile(const std::string &fn,
                      const IO::Plugins &p=IO::defaultPlugins());
    // read with explicit format
    IO::Data readFile(const std::string &fn, const IO::Plugin* plug);
    bool     writeFile(const std::string &fn, const IO::Plugin* plug, const Molecule &m,
                       size_t idx=-1ul,
                       const IO::BaseParam *p=nullptr,
                       const IO::BasePreset *c=nullptr);
    std::optional<const IO::Plugin*> guessFmt(std::string fn,
                                   IO::Plugin::Args arg = IO::Plugin::Read,
                                   const IO::Plugins &p=IO::defaultPlugins());
}

#endif // IOWRAPPER

