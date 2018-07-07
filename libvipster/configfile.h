#ifndef CONFIG
#define CONFIG

#include <string>
#include <map>
#include <array>
//#include <variant>
#include <memory>
#include <cstdlib>
#ifdef __APPLE__
#include <CoreFoundation/CoreFoundation.h>
#endif
#ifdef _WIN32
#include <windows.h>
#endif
#ifdef __linux__
#include <unistd.h>
#include <libgen.h>
#endif

namespace Vipster{

// call this function to initialize the library-state
bool readConfig();
bool saveConfig();

#ifdef __EMSCRIPTEN__
const std::string sys_path = "/";
const std::string sys_config = "/vipster.json";
const std::string user_path{}; // not used
const std::string user_config{}; // not used
#elif __linux__
#ifdef RELPATH
const std::string sys_path = [](){
    char path[4096]; // let's hope this is big enough
    readlink("/proc/self/exe", path, 4096);
    return std::string{dirname(dirname(path))}+"/share/vipster";
}();
#else
#ifndef PREFIX
#define PREFIX /usr
#endif
#define TO_STR2(x) #x
#define TO_STR(x) TO_STR2(x)
const std::string sys_path = std::string{TO_STR(PREFIX)} + "/share/vipster";
#endif
const std::string sys_config = sys_path + "/default.json";
const std::string user_path = std::getenv("HOME");
const std::string user_config = user_path + "/.vipster.json";
#elif _WIN32
const std::string sys_path = [](){
    const HMODULE hModule = GetModuleHandleA(nullptr);
    CHAR path[MAX_PATH];
    GetModuleFileNameA(hModule, path, MAX_PATH);
    std::string temp{path};
    auto pos = temp.find_last_of('\\');
    return temp.substr(0, pos);
}();
const std::string sys_config = sys_path + "\\default.json";
const std::string user_path = std::string{std::getenv("APPDATA")} + "\\vipster";
const std::string user_config = user_path + "\\vipster.json";
#elif __APPLE__
const std::string sys_path = [](){
    CFURLRef appUrlRef = CFBundleCopyBundleURL(CFBundleGetMainBundle());
    CFStringRef macPath = CFURLCopyFileSystemPath(appUrlRef, kCFURLPOSIXPathStyle);
    const char *pathPtr = CFStringGetCStringPtr(macPath, CFStringGetSystemEncoding());
    CFRelease( appUrlRef );
    CFRelease( macPath );
    return std::string{pathPtr};
}();
const std::string sys_config = sys_path + "/Contents/Resources/default.json";
const std::string user_path = std::string{std::getenv("HOME")} + "/Library/Application Support/vipster";
const std::string user_config = user_path + "/vipster.json";
#endif

}

#endif // CONFIG
