#ifndef CONFIG
#define CONFIG

#include <string>
#include <map>
#include <array>
#ifdef _WIN32
#include <windows.h>
#endif

namespace Vipster{

#ifdef __EMSCRIPTEN__
const std::string user_config = "/vipster.json";
#elif __linux__
#ifndef PREFIX
#define PREFIX /usr/share/
#endif
const std::string sys_config = "PREFIXvipster.json";
const std::string user_config = std::string(std::getenv("HOME"))+"/.vipster.json";
#elif _WIN32
const std::string sys_config = [](){
    const HMODULE hModule = GetModuleHandleA(nullptr);
    CHAR path[MAX_PATH];
    GetModuleFileNameA(hModule, path, MAX_PATH);
    std::string temp{path};
    auto pos = temp.find_last_of('\\');
    return temp.substr(0, pos).append("\\default.json");
}();
const std::string user_config = std::string(std::getenv("USERPROFILE"))+"/vipster.json";
#elif __APPLE__
//TODO: should be something like ~/Library/Application Support/vipster/vipster.json
//TODO: sys_config maybe Library/Application Support/vipster/vipster.json
const std::string user_config = std::string(std::getenv("HOME"))+"/.vipster.json";
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

PseMap readPse(void);

extern const PseMap pse;

}

#endif // CONFIG

