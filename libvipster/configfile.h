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

namespace Vipster{

// call this functions to read user-defined settings
bool readConfig();
bool saveConfig();

}

#endif // CONFIG
