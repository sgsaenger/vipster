#include "configfile.h"
#include "io.h"
#include "json.hpp"

#include <fstream>
#include <tuple>

using json = nlohmann::json;
using namespace Vipster;

std::tuple<PeriodicTable, Settings, IO::Parameters, IO::Configs> Vipster::readConfig()
{
    ConfigState retVal{};
    std::ifstream conf_file{user_config};
    PeriodicTable &pte = std::get<0>(retVal);
    Settings &settings = std::get<1>(retVal);
    IO::Parameters &params = std::get<2>(retVal);
    IO::Configs &configs = std::get<3>(retVal);
    // load defaults as a minimum starting point
    pte = Vipster::pte;
    settings = Vipster::settings;
    params = IO::defaultParams();
    configs = IO::defaultConfigs();
    if(conf_file){
        try {
            json loc_file;
            conf_file >> loc_file;
            // Periodic table
            if(loc_file.find("PSE") != loc_file.end()){
                from_json(loc_file["PSE"], pte);
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
                            auto& tmp = params[pair.first];
                            for(auto it2=it.value().begin(); it2!=it.value().end(); ++it2){
                                tmp[it2.key()] = plugin->makeParam(it2.key());
                                tmp[it2.key()]->parseJson(it2);
                            }
                            continue;
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
                            auto& tmp = configs[pair.first];
                            for(auto it2=it.value().begin(); it2!=it.value().end(); ++it2){
                                tmp[it2.key()] = plugin->makeConfig(it2.key());
                                tmp[it2.key()]->parseJson(it2);
                            }
                            continue;
                        }
                    }
                }
            }
        } catch (const json::exception& e) {
            std::cout << e.what() << std::endl;
        }
    }
    return retVal;
}

void Vipster::saveConfig(const ConfigState& cs)
{
    std::ofstream conf_file{user_config};
    if(!conf_file){
        throw IO::Error("Cannot open config file at \""+user_config+"\" for writing");
    }
    const PeriodicTable &pte = std::get<0>(cs);
    const Settings &settings = std::get<1>(cs);
    const IO::Parameters &params = std::get<2>(cs);
    const IO::Configs &configs = std::get<3>(cs);
    json json_buf;
    json_buf["PSE"] = pte;
    json_buf["Settings"] = settings;
    json_buf["Parameters"] = json{};
    for(const auto& plug: IOPlugins){
        if(plug.second->arguments & IO::Plugin::Args::Param){
            json& j = json_buf["Parameters"][plug.second->command];
            for(const auto& p: params.at(plug.first)){
                j[p.first] = p.second->toJson();
            }
        }
    }
    json_buf["Configs"] = json{};
    for(const auto& plug: IOPlugins){
        if(plug.second->arguments & IO::Plugin::Args::Config){
            json& j = json_buf["Configs"][plug.second->command];
            for(const auto& p: configs.at(plug.first)){
                j[p.first] = p.second->toJson();
            }
        }
    }
    conf_file << json_buf.dump(2);
    return;
}
