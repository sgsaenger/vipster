#ifndef CONFIG
#define CONFIG

#include "json.hpp"

#include <string>
#include <map>
#include <array>
#include <variant>
#include <memory>
#include <cstdlib>
#ifdef _WIN32
#include <windows.h>
#endif

namespace Vipster{

#ifdef __EMSCRIPTEN__
const std::string sys_path = "/";
const std::string sys_config = "/vipster.json";
const std::string user_path{};
const std::string user_config{};
#elif __linux__
#ifndef PREFIX
#define PREFIX /usr/share/
#endif
const std::string sys_path = "PREFIXvipster";
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
//TODO: maybe Library/Application Support/vipster/vipster.json? or in bundle?
const std::string sys_path = "";
const std::string sys_config = "";
//TODO: should be something like ~/Library/Application Support/vipster/vipster.json
const std::string user_path = std::string{std::getenv("HOME")} + "/Library/Application Support/vipster";
const std::string user_config = user_path + "/.vipster.json";
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

struct BaseParam
{
    std::string name;
    virtual std::unique_ptr<BaseParam> copy() = 0;
    virtual ~BaseParam() = default;
};

struct BaseConfig
{
    std::string name;
    virtual std::unique_ptr<BaseConfig> copy() = 0;
    virtual ~BaseConfig() = default;
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

using Settings = std::map<std::string, std::variant<bool, int, float, std::string>>;
using Parameters = std::multimap<IOFmt, std::unique_ptr<BaseParam>>;

extern PseMap pse;
extern Settings settings;
extern Parameters params;

}

#endif // CONFIG
