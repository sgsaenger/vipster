#ifndef CONFIG
#define CONFIG

#include "json.hpp"
#include "global.h"
#include "bond.h"

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

using ColVec = std::array<uint8_t, 4>;

struct PseEntry{
    std::string     PWPP;
    std::string     CPPP;
    std::string     CPNL;
    unsigned int    Z;
    float           m;
    float           bondcut;
    float           covr;
    float           vdwr;
    ColVec          col;
};

enum class IOFmt{XYZ, PWI, PWO, LMP, DMP};

class BaseParam
{
public:
    std::string name;
    virtual std::unique_ptr<BaseParam> copy() = 0;
    virtual ~BaseParam() = default;
protected:
    BaseParam(std::string);
};

class BaseConfig
{
public:
    std::string name;
    virtual std::unique_ptr<BaseConfig> copy() = 0;
    virtual ~BaseConfig() = default;
protected:
    BaseConfig(std::string);
};


class PseMap:private std::map<std::string,PseEntry>
{
public:
    using std::map<std::string,PseEntry>::begin;
    using std::map<std::string,PseEntry>::end;
    using std::map<std::string,PseEntry>::at;
    using std::map<std::string,PseEntry>::emplace;
    using std::map<std::string,PseEntry>::find;
    using std::map<std::string,PseEntry>::size;
    PseMap(bool r=false):root(r){}
    PseEntry& operator [](const std::string &k);
private:
    bool root;
};

template<typename T>
struct Setting{
    using ValueType = T;
    std::string name;
    T val;
};
struct Settings{
    Setting<bool> atRadVdW{"Atom radius VdW", false};
    Setting<float> atRadFac{"Atom radius factor", bohrrad};
    Setting<float> bondRad{"Bond radius", bohrrad};
    Setting<float> bondCutFac{"Bond cutoff factor", 1.1f};
    Setting<BondFrequency> bondFreq{"Bond frequency", BondFrequency::Always};
    Setting<BondLevel> bondLvl{"Bond level", BondLevel::Cell};
    Setting<bool> showBonds{"Show bonds", true};
    Setting<bool> showCell{"Show cell", true};
    Setting<bool> antialias{"Antialiasing", true};
    Setting<bool> perspective{"Perspective projection", false};
    Setting<ColVec> selCol{"Selection color", ColVec{0, 0, 80, 80}};
    Setting<std::string> PWPP{"Default PWScf PP-suffix", ""};
    Setting<std::string> CPPP{"Default CPMD PP-suffix", ""};
    Setting<std::string> CPNL{"Default CPMD Nonlocality", "LMAX=F"};
};
//TODO: void to_json(nlohmann::json& j,const Settings& s);
void from_json(const nlohmann::json& j, Settings& s);

using Parameters = std::multimap<IOFmt, std::unique_ptr<BaseParam>>;
using Configs = std::multimap<IOFmt, std::unique_ptr<BaseConfig>>;

extern PseMap pse;
extern Settings settings;
extern Parameters params;
extern Configs configs;

}

#endif // CONFIG
