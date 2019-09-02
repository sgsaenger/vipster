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
#include "io/parameters.h"
#include "io/configs.h"

namespace Vipster{

// call this functions to read user-defined settings
using ConfigState = std::tuple<PeriodicTable, Settings, IO::Parameters, IO::Configs>;
ConfigState readConfig();
void saveConfig(const ConfigState &);

#if __linux__ || defined(__FreeBSD__)
const std::string user_path = [](){
        auto tmp = std::getenv("XDG_CONFIG_HOME");
        return tmp == nullptr ? std::string{std::getenv("HOME")}+"/.config" : tmp;
    }();
const std::string user_config = user_path + "/vipster.json";
#elif _WIN32
const std::string user_path = std::string{std::getenv("APPDATA")} + "\\vipster";
const std::string user_config = user_path + "\\vipster.json";
#elif __APPLE__
const std::string user_path = std::string{std::getenv("HOME")} + "/Library/Application Support/vipster";
const std::string user_config = user_path + "/vipster.json";
#else
const std::string user_path{}; // not used
const std::string user_config{}; // not used
#endif

}

#endif // CONFIG
