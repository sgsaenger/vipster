#ifndef CONFIG
#define CONFIG

#include <tuple>
#include <filesystem>

#include "periodictable.h"
#include "settings.h"
#include "plugin.h"
#include "parameters.h"
#include "presets.h"

namespace Vipster{

// call this functions to read user-defined settings
using ConfigState = std::tuple<PeriodicTable,
                               Settings,
                               PluginList,
                               ParameterMap,
                               PresetMap
                               >;
ConfigState readConfig();
void saveConfig(const ConfigState &);
std::filesystem::path getConfigDir();

}

#endif // CONFIG
