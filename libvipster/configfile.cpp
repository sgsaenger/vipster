#include "configfile.h"
#include "iowrapper.h"
#include "pse.h"
#include "settings.h"
#include "parameters.h"
#include "configs.h"
#include "json.hpp"

#include <fstream>
#include <tuple>

using json = nlohmann::json;
using namespace Vipster;

bool readConfig();

namespace Vipster{
PseMap pse;
Settings settings;
Parameters params;
Configs configs;
static bool config_loaded = readConfig();
}

bool readConfig()
{
    PseMap pse{true};
    Settings settings;
    Parameters params;
    Configs configs;
    std::ifstream conf_file{user_config};
    if(!conf_file){
        conf_file = std::ifstream{sys_config};
    }
    if(conf_file){
        json loc_file;
        conf_file >> loc_file;
        // PSE
        for(auto it=loc_file["PSE"].begin();it!=loc_file["PSE"].end();++it)
        {
            auto v = it.value();
            pse.emplace(it.key(),PseEntry{v["PWPP"],
                    v["CPPP"],v["CPNL"],v["Z"],
                    v["m"],v["bondcut"],v["covr"],
                    v["vdwr"],v["col"]});
        }
        // General settings
        settings = loc_file["Settings"];
        // Parameter sets
        for(const auto& pw:loc_file["Parameters"]["PWI"]){
            params.emplace(IOFmt::PWI, std::make_unique<IO::PWParam>(pw.get<IO::PWParam>()));
        }
    }
    // ensure fallback-values are present
    if(pse.find("")==pse.end()){
        pse.emplace("", PseEntry{"","","",0,0,0,1.46f,3.21f,{{0,0,0,255}}});
    }
    if(configs.find(IOFmt::XYZ) == configs.end()){
        configs.emplace(IOFmt::XYZ, std::make_unique<IO::XYZConfig>(IO::XYZConfigDefault));
    }
    if(params.find(IOFmt::PWI) == params.end()){
        params.emplace(IOFmt::PWI, std::make_unique<IO::PWParam>(IO::PWParamDefault));
    }
    if(configs.find(IOFmt::PWI) == configs.end()){
        configs.emplace(IOFmt::PWI, std::make_unique<IO::PWConfig>(IO::PWConfigDefault));
    }
    if(configs.find(IOFmt::LMP) == configs.end()){
        configs.emplace(IOFmt::LMP, std::make_unique<IO::LmpConfig>(IO::LmpConfigDefault));
    }
    Vipster::pse = std::move(pse);
    Vipster::settings = std::move(settings);
    Vipster::params = std::move(params);
    Vipster::configs = std::move(configs);
    return true;
}
