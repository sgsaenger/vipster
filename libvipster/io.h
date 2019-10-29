#ifndef IOWRAPPER
#define IOWRAPPER

#include <map>

#include "molecule.h"
#include "io/data.h"
#include "io/plugin.h"
#include "io/parameters.h"
#include "io/configs.h"

//TODO: check std::ios_base::sync_with_stdio(false)
namespace Vipster{
    // convenience containers and defaults
    namespace IO{
        using Plugins = std::vector<const Plugin*>;
        Plugins defaultPlugins();
        using Configs = std::map<const Plugin*, std::map<std::string, std::unique_ptr<BaseConfig>>>;
        Configs defaultConfigs(const Plugins& p);
        using Parameters = std::map<const Plugin*, std::map<std::string, std::unique_ptr<BaseParam>>>;
        Parameters defaultParams(const Plugins& p);
    }

    // read with format guess
    IO::Data readFile(const std::string &fn,
                      const IO::Plugins &p=IO::defaultPlugins());
    // read with explicit format
    IO::Data readFile(const std::string &fn, const IO::Plugin* plug);
    bool     writeFile(const std::string &fn, const IO::Plugin* plug, const Molecule &m,
                       const IO::BaseParam *p=nullptr,
                       const IO::BaseConfig *c=nullptr,
                       size_t idx=-1ul);
    std::optional<const IO::Plugin*> guessFmt(std::string fn,
                                   IO::Plugin::Args arg = IO::Plugin::Read,
                                   const IO::Plugins &p=IO::defaultPlugins());
}

#endif // IOWRAPPER

