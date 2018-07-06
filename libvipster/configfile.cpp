#include "configfile.h"
#include "io.h"
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
}

bool Vipster::readConfig()
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
        for(auto it=loc_file["Parameters"].begin(); it!=loc_file["Parameters"].end();++it)
        {
            for(const auto& pair: IOPlugins){
                const auto& plugin = pair.second;
                if(it.key() == plugin->command){
                    for(const auto& entry: it.value()){
                        auto tmp = plugin->makeParam("");
                        tmp->parseJson(entry);
                        params.emplace(pair.first, std::move(tmp));
                    }
                }
            }
        }
    }
    // ensure fallback-values are present
    if(pse.find("")==pse.end()){
        pse.emplace("", PseEntry{"","","",0,0,0,1.46f,3.21f,{{0,0,0,255}}});
    }
    for(const auto& pair: IOPlugins){
        auto fmt = pair.first;
        const auto& plugin = pair.second;
        if(plugin->arguments & IO::Plugin::Args::Config &&
           configs.find(fmt) == configs.end()){
           configs.emplace(fmt, plugin->makeConfig("default"));
        }
        if(plugin->arguments & IO::Plugin::Args::Param &&
           params.find(fmt) == params.end()){
           params.emplace(fmt, plugin->makeParam("default"));
        }
    }
    Vipster::pse = std::move(pse);
    Vipster::settings = std::move(settings);
    Vipster::params = std::move(params);
    Vipster::configs = std::move(configs);
    return true;
}
