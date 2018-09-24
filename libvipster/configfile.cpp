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

bool Vipster::readConfig()
{
    std::ifstream conf_file{user_config};
    if(conf_file){
        json loc_file;
        conf_file >> loc_file;
        // PSE
        if(loc_file.find("PSE") != loc_file.end()){
            from_json(loc_file["PSE"], pse);
        }
        // General settings
        if(loc_file.find("Settings") != loc_file.end()){
            from_json(loc_file["Settings"], settings);
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
