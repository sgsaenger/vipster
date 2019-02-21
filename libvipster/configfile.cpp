#include "configfile.h"
#include "io.h"
#include "periodictable.h"
#include "settings.h"
#include "io/parameters.h"
#include "io/presets.h"
#include "json.hpp"

#include <fstream>
#include <tuple>
#include <filesystem>

using json = nlohmann::json;
using namespace Vipster;
namespace fs = std::filesystem;

fs::path getConfigDir(){
#if defined(__linux__) || defined(__FreeBSD__)
    auto tmp = std::getenv("XDG_CONFIG_HOME");
    return fs::path{tmp == nullptr ? std::string{std::getenv("HOME")}+"/.config" : tmp}/"vipster";
#elif _WIN32
    return fs::path{std::getenv("APPDATA")}/"vipster";
#elif __APPLE__
    return fs::path{std::getenv("HOME")}/"Library/Application Support/vipster";
#endif
}

bool Vipster::readConfig()
{
    bool success{true};
    auto dir = getConfigDir();
    if(!fs::exists(dir)){
        std::cout << "Config directory at \"" << dir << "\" does not exist" << std::endl;
        return false;
    }
    // PSE
    fs::path pse_path = dir/"pse.json";
    std::ifstream pse_file{pse_path};
    if(pse_file){
        try {
            json j;
            pse_file >> j;
            pse = j;
        } catch (const json::exception& e) {
            std::cout << "Error when reading PSE: " << e.what() << std::endl;
            success = false;
        }
    }
    // Settings
    fs::path settings_path = dir/"settings.json";
    std::ifstream settings_file{settings_path};
    if(settings_file){
        try {
            json j;
            settings_file >> j;
            settings = j;
        } catch (const json::exception& e) {
            std::cout << "Error when reading settings: " << e.what() << std::endl;
            success = false;
        }
    }
    // Parameter sets
    fs::path param_path = dir/"parameters.json";
    std::ifstream param_file{param_path};
    if(param_file){
        try {
            json j;
            param_file >> j;
            for(const auto& pair: IOPlugins){
                const auto& plugin = pair.second;
                if(!(plugin->arguments & IO::Plugin::Param)) continue;
                auto pos = j.find(plugin->command);
                if(pos != j.end()){
                    auto& tmp = params[pair.first];
                    for(auto param: pos->items()){
                        tmp[param.key()] = plugin->makeParam();
                        tmp[param.key()]->parseJson(param.value());
                    }
                }
            }
        } catch (const json::exception& e) {
            std::cout << "Error when reading parameters: " << e.what() << std::endl;
            success = false;
        }
    }
    // IO-presets
    fs::path presets_path = dir/"iopresets.json";
    std::ifstream presets_file{presets_path};
    if(presets_file){
        try {
            json j;
            presets_file >> j;
            for(const auto& pair: IOPlugins){
                const auto& plugin = pair.second;
                if(!(plugin->arguments & IO::Plugin::Preset)) continue;
                auto pos = j.find(plugin->command);
                if(pos != j.end()){
                    auto& tmp = presets[pair.first];
                    for(auto preset: pos->items()){
                        tmp[preset.key()] = plugin->makePreset();
                        tmp[preset.key()]->parseJson(preset.value());
                    }
                }
            }
        } catch (const json::exception& e) {
            std::cout << "Error when reading IO-presets: " << e.what() << std::endl;
            success = false;
        }
    }
    return success;
}

bool Vipster::saveConfig()
{
    auto dir = getConfigDir();
    bool success{true};
    json j;
    if(!fs::exists(dir)){
        std::error_code error;
        fs::create_directories(dir, error);
        if(error){
            std::cout << "Couldn't create directory \"" << dir
                      << "\" to write configuration" << std::endl;
            return false;
        }
    }
    // PSE
    fs::path pse_path = dir/"pse.json";
    std::ofstream pse_file{pse_path};
    if(!pse_file){
        std::cout << "Can not open file at \"" << std::string{pse_path} << "\" for writing PSE" << std::endl;
        success = false;
    }
    j = pse;
    pse_file << j.dump(2);
    // Settings
    fs::path settings_path = dir/"settings.json";
    std::ofstream settings_file{settings_path};
    if(!pse_file){
        std::cout << "Can not open file at \"" << std::string{settings_path} << "\" for writing settings";
        success = false;
    }
    j = settings;
    settings_file << j.dump(2);
    // Parameters
    fs::path param_path = dir/"parameters.json";
    std::ofstream param_file{param_path};
    if(!pse_file){
        std::cout << "Can not open file at \"" << std::string{param_path} << "\" for writing parameter sets";
        success = false;
    }
    j = json{};
    for(const auto& pair: params){
        const auto& plugin = IOPlugins.at(pair.first);
        if(!(plugin->arguments & IO::Plugin::Param)) continue;
        const auto& com = plugin->command;
        j[com] = json{};
        auto& tmp = j[com];
        for(const auto& param: pair.second){
            tmp[param.first] = *param.second;
        }
    }
    param_file << j.dump(2);
    // IO-Presets
    fs::path preset_path = dir/"iopresets.json";
    std::ofstream preset_file{preset_path};
    if(!pse_file){
        std::cout << "Can not open file at \"" << std::string{preset_path} << "\" for writing IO-presets";
        success = false;
    }
    j = json{};
    for(const auto& pair: presets){
        const auto& plugin = IOPlugins.at(pair.first);
        if(!(plugin->arguments & IO::Plugin::Preset)) continue;
        const auto& com = plugin->command;
        j[com] = json{};
        auto& tmp = j[com];
        for(const auto& preset: pair.second){
            tmp[preset.first] = *preset.second;
        }
    }
    preset_file << j.dump(2);
    return success;
}
