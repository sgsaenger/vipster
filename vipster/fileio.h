#ifndef IOWRAPPER
#define IOWRAPPER

#include <map>
#include <filesystem>
#include <optional>

#include "molecule.h"
#include "iotuple.h"
#include "plugin.h"
#include "parameters.h"
#include "presets.h"

//TODO: check std::ios_base::sync_with_stdio(false)
namespace Vipster{
    // read with format guess
    IOTuple readFile(const std::string &fn,
                      const PluginList &p=defaultPlugins());
    // read with explicit format
    IOTuple readFile(const std::string &fn, const Plugin* plug);
    bool    writeFile(const std::string &fn, const Plugin* plug, const Molecule &m,
                      std::optional<size_t> idx={},
                      const std::optional<Parameter>& p=std::nullopt,
                      const std::optional<Preset>& c=std::nullopt);
    const Plugin* guessFmt(std::string fn,
                           const PluginList &p=defaultPlugins());
    // RAII wrapper for temp folder
    namespace detail {
        class TempWrap{
        public:
            TempWrap();
            const std::filesystem::path& getPath() const;
        private:
            TempWrap(const TempWrap&) = delete;
            std::filesystem::path tmppath;
        };
        extern const TempWrap tempwrap;
    }
    const std::filesystem::path& getTempPath();
}

#endif // IOWRAPPER

