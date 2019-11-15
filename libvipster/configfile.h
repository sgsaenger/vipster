#ifndef CONFIG
#define CONFIG

#include <tuple>

#include "periodictable.h"
#include "settings.h"
#include "fileio.h"

namespace Vipster{

// call this functions to read user-defined settings
using ConfigState = std::tuple<PeriodicTable,
                               Settings,
                               IO::Plugins,
                               IO::Parameters,
                               IO::Presets
                               >;
ConfigState readConfig();
void saveConfig(const ConfigState &);

}

#endif // CONFIG
