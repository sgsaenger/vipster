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

/* obtain platform-specific path for config
 * Windows: %APPDATA%/vipster
 * macOS: $HOME/Library/Application Support/vipster
 * Linux/other: $XDG_CONFIG_HOME/vipster or $HOME/.config/vipster
 */
std::filesystem::path getConfigDir();

// Read user-specified settings or defaults (call this on app-start)
struct ConfigState{
    PeriodicTable   periodicTable;
    Settings        settings;
    PluginList      plugins;
    ParameterMap    parameters;
    PresetMap       presets;
};
ConfigState readConfig();

// Store user-specified settings (call this on app-shutdown)
void saveConfig(const ConfigState &);

}

#endif // CONFIG
