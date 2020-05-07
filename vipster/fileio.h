#ifndef IOWRAPPER
#define IOWRAPPER

#include <map>
#include <filesystem>
#include <optional>

#include "molecule.h"
#include "io/data.h"
#include "io/plugin.h"
#include "io/parameters.h"
#include "io/presets.h"

//TODO: check std::ios_base::sync_with_stdio(false)
namespace Vipster{
    // read with format guess
    IO::Data readFile(const std::string &fn,
                      const IO::Plugins &p=IO::defaultPlugins());
    // read with explicit format
    IO::Data readFile(const std::string &fn, const IO::Plugin* plug);
    bool     writeFile(const std::string &fn, const IO::Plugin* plug, const Molecule &m,
                       std::optional<size_t> idx={},
                       const std::optional<IO::Parameter>& p=std::nullopt,
                       const std::optional<IO::Preset>& c=std::nullopt);
    const IO::Plugin* guessFmt(std::string fn,
                               const IO::Plugins &p=IO::defaultPlugins());
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

