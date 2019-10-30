#ifndef CONFIG
#define CONFIG

#include <string>
#include <map>
#include <array>
#include <memory>
#include <cstdlib>
#ifdef __APPLE__
#include <CoreFoundation/CoreFoundation.h>
#endif
#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#endif
#ifdef __linux__
#include <unistd.h>
#include <libgen.h>
#endif

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
