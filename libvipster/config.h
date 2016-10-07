#ifndef CONFIG
#define CONFIG

#include "definitions.h"

namespace Vipster{

#ifdef __linux__
#ifndef PREFIX
#define PREFIX /usr/share/
#endif
const std::string sys_config = "PREFIXvipster.json";
const std::string user_config = std::string(std::getenv("HOME"))+"/.vipster.json";
#elif _WIN32
//WIP
HMODULE hModule = GetModuleHandleW(NULL);
WCHAR path[MAX_PATH];
GetModuleFileNameW(hModule,path,MAX_PATH);
const std::string sys_config = "./vipster.json";
const std::string user_config = std::string(std::getenv("USERPROFILE"))+"/vipster.json";
#endif

struct PseEntry{
    std::string PWPP;
    std::string CPPP;
    std::string CPNL;
    uint        Z;
    float       m;
    float       bondcut;
    float       covr;
    float       vdwr;
    std::vector<float> col;
};

class PseMap:public std::map<std::string,PseEntry>
{
public:
    PseMap(const std::map<std::string,PseEntry> *x):internal(x){};
    PseEntry& operator [](const std::string &k);
private:
    const std::map<std::string,PseEntry> * internal;
};

std::map<std::string,PseEntry> readPse(void);

const std::map<std::string,PseEntry> pse = readPse();

}

#endif // CONFIG

