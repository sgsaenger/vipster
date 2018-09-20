#include "configfile.h"
#include "io.h"
#include "pse.h"
#include "settings.h"
#include "io/parameters.h"
#include "io/configs.h"
#include "json.hpp"

#include <fstream>
#include <tuple>

using json = nlohmann::json;
using namespace Vipster;

bool readConfig();

namespace Vipster{
PseMap pse;
Settings settings;
IO::Parameters params;
IO::Configs configs;
}

bool Vipster::readConfig()
{
    PseMap pse{true};
    Settings settings;
    IO::Parameters params;
    IO::Configs configs;
    std::ifstream conf_file{user_config};
    if(!conf_file){
        conf_file = std::ifstream{sys_config};
    }
    if(conf_file){
        json loc_file;
        conf_file >> loc_file;
        // PSE
        if(loc_file.find("PSE") != loc_file.end()){
            pse = loc_file["PSE"];
        }
        // General settings
        if(loc_file.find("Settings") != loc_file.end()){
            settings = loc_file["Settings"];
        }
        // Parameter sets
        if(loc_file.find("Parameters") != loc_file.end()){
            for(auto it=loc_file["Parameters"].begin(); it!=loc_file["Parameters"].end(); ++it)
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
        // Config presets
        if(loc_file.find("Configs") != loc_file.end()){
            for(auto it=loc_file["Configs"].begin(); it!=loc_file["Configs"].end(); ++it){
                for(const auto& pair: IOPlugins){
                    const auto& plugin = pair.second;
                    if(it.key() == plugin->command){
                        for(const auto& entry: it.value()){
                            auto tmp = plugin->makeConfig("");
                            tmp->parseJson(entry);
                            configs.emplace(pair.first, std::move(tmp));
                        }
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

bool Vipster::saveConfig()
{
    std::ofstream conf_file{user_config};
    if(!conf_file){
        throw IO::Error("Cannot open config file at \""+user_config+"\" for writing");
    }
    json json_buf;
    json_buf["PSE"] = pse;
    json_buf["Settings"] = settings;
    json_buf["Parameters"] = json{};
    for(const auto& plug: IOPlugins){
        if(plug.second->arguments & IO::Plugin::Args::Param){
            json_buf["Parameters"][plug.second->command] = json::array();
            json& j = json_buf["Parameters"][plug.second->command];
            for(const auto& p: params){
                if(p.first == plug.first){
                    j.push_back(p.second->toJson());
                }
            }
        }
    }
    json_buf["Configs"] = json{};
    for(const auto& plug: IOPlugins){
        if(plug.second->arguments & IO::Plugin::Args::Config){
            json_buf["Configs"][plug.second->command] = json::array();
            json& j = json_buf["Configs"][plug.second->command];
            for(const auto& p: configs){
                if(p.first == plug.first){
                    j.push_back(p.second->toJson());
                }
            }
        }
    }
    conf_file << json_buf.dump(2);
    return true;
}
