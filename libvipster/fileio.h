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
        using Presets = std::map<const Plugin*, std::map<std::string, BasePreset>>;
        Presets defaultPresets(const Plugins& p);
        using Parameters = std::map<const Plugin*, std::map<std::string, BaseParam>>;
        Parameters defaultParams(const Plugins& p);
    }

    // read with format guess
    IO::Data readFile(const std::string &fn,
                      const IO::Plugins &p=IO::defaultPlugins());
    // read with explicit format
    IO::Data readFile(const std::string &fn, const IO::Plugin* plug);
    bool     writeFile(const std::string &fn, const IO::Plugin* plug, const Molecule &m,
                       std::optional<size_t> idx={},
                       const std::optional<IO::BaseParam>& p=std::nullopt,
                       const std::optional<IO::BasePreset>& c=std::nullopt);
    std::optional<const IO::Plugin*> guessFmt(std::string fn,
                                   const IO::Plugins &p=IO::defaultPlugins());
}

#endif // IOWRAPPER

