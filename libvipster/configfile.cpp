#include "configfile.h"
#include "nlohmann/json.hpp"

#include <fstream>

#include <cstdlib>
#if defined(__APPLE__)
#include <CoreFoundation/CoreFoundation.h>
#include <dlfcn.h>
#elif defined(_WIN32)
#define NOMINMAX
#include <windows.h>
#else
#include <unistd.h>
#include <libgen.h>
#include <dlfcn.h>
#endif

using json = nlohmann::json;
using namespace Vipster;
namespace fs = std::filesystem;

#if defined(_WIN32)
#define LIBEXTENSION ".dll"
#elif defined(__APPLE__)
#define LIBEXTENSION ".dylib"
#else
#define LIBEXTENSION ".so"
#endif

const IO::Plugin* openPlugin(fs::path name)
{
#if defined(_WIN32)
    try{
        auto file = LoadLibrary(std::string{name.string()}.c_str());
        if(file == NULL) return nullptr;
        return reinterpret_cast<const IO::Plugin*>(GetProcAddress(file, "plugin"));
    }catch(...){
        return nullptr;
    }
#else
    try{
        auto *file = dlopen(name.c_str(), RTLD_LAZY);
        if(!file) return nullptr;
        return static_cast<const IO::Plugin*>(dlsym(file, "plugin"));
    }catch(...){
        return nullptr;
    }
#endif
}

fs::path Vipster::getConfigDir(){
#if defined(__linux__) || defined(__FreeBSD__)
    auto tmp = std::getenv("XDG_CONFIG_HOME");
    return fs::path{tmp == nullptr ? std::string{std::getenv("HOME")}+"/.config" : tmp}/"vipster";
#elif defined(_WIN32)
    return fs::path{std::getenv("APPDATA")}/"vipster";
#elif defined(__APPLE__)
    return fs::path{std::getenv("HOME")}/"Library/Application Support/vipster";
#endif
}

ConfigState Vipster::readConfig()
{
    // create state
    ConfigState retVal{};
    PeriodicTable &pte = std::get<0>(retVal);
    Settings &settings = std::get<1>(retVal);
    IO::Plugins &plugins = std::get<2>(retVal);
    IO::Parameters &params = std::get<3>(retVal);
    IO::Presets &presets = std::get<4>(retVal);
    // load defaults as a minimum starting point
    pte = Vipster::pte;
    settings = Vipster::settings;
    plugins = IO::defaultPlugins();
    params = IO::defaultParams(plugins);
    presets = IO::defaultPresets(plugins);
    // make sure config dir is available
    auto dir = getConfigDir();
    if(!fs::exists(dir)){
        std::cerr << "Config directory at \"" << dir << "\" does not exist" << std::endl;
        return retVal;
    }
    // User-created plugins
    auto pluginDir = dir/"plugins";
    if(fs::exists(pluginDir)){
        for(const auto& file: fs::directory_iterator(pluginDir)){
            if(file.path().extension() != LIBEXTENSION) continue;
            auto* plug = openPlugin(file);
            if(plug){
                std::cerr << "Loading plugin " << file << std::endl;
                plugins.push_back(plug);
            }
        }
    }
    // Periodic table
    fs::path pte_path = dir/"periodictable.json";
    std::ifstream pte_file{pte_path};
    if(pte_file){
        try {
            json j;
            pte_file >> j;
            pte = j;
        } catch (const json::exception& e) {
            std::cerr << "Error when reading Periodic table: " << e.what() << std::endl;
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
            std::cerr << "Error when reading settings: " << e.what() << std::endl;
        }
    }
    // Parameter sets
    fs::path param_path = dir/"parameters.json";
    std::ifstream param_file{param_path};
    if(param_file){
        try {
            json j;
            param_file >> j;
            for(const auto& plugin: plugins){
                if(!plugin->makeParam) continue;
                auto pos = j.find(plugin->command);
                if(pos != j.end()){
                    auto& tmp = params[plugin];
                    for(auto param: pos->items()){
                        tmp[param.key()] = plugin->makeParam();
                        tmp[param.key()]->parseJson(param.value());
                    }
                }
            }
        } catch (const json::exception& e) {
            std::cerr << "Error when reading parameters: " << e.what() << std::endl;
        }
    }
    // IO-presets
    fs::path presets_path = dir/"iopresets.json";
    std::ifstream presets_file{presets_path};
    if(presets_file){
        try {
            json j;
            presets_file >> j;
            for(const auto& plugin: plugins){
                if(!plugin->makePreset) continue;
                auto pos = j.find(plugin->command);
                if(pos != j.end()){
                    auto& tmp = presets[plugin];
                    for(auto preset: pos->items()){
                        tmp[preset.key()] = plugin->makePreset();
                        tmp[preset.key()]->parseJson(preset.value());
                    }
                }
            }
        } catch (const json::exception& e) {
            std::cerr << "Error when reading IO-presets: " << e.what() << std::endl;
        }
    }
    return retVal;
}

void Vipster::saveConfig(const ConfigState& cs)
{
    // convenience refs
    const PeriodicTable &pte = std::get<0>(cs);
    const Settings &settings = std::get<1>(cs);
    const IO::Plugins &plugins = std::get<2>(cs);
    const IO::Parameters &params = std::get<3>(cs);
    const IO::Presets &presets = std::get<4>(cs);
    // check dir, try to create if it doesn't exist
    auto dir = getConfigDir();
    bool success{true};
    json j;
    if(!fs::exists(dir)){
        std::error_code error;
        fs::create_directories(dir, error);
        if(error){
            std::cerr << "Couldn't create directory \"" << dir
                      << "\" to write configuration" << std::endl;
            return;
        }
    }
    // Periodic table
    fs::path pte_path = dir/"periodictable.json";
    std::ofstream pte_file{pte_path};
    if(!pte_file){
        std::cerr << "Can not open file at \"" << pte_path << "\" for writing Periodic Table" << std::endl;
        success = false;
    }
    j = pte;
    pte_file << j.dump(2);
    // Settings
    fs::path settings_path = dir/"settings.json";
    std::ofstream settings_file{settings_path};
    if(!settings_file){
        std::cerr << "Can not open file at \"" << settings_path << "\" for writing settings";
        success = false;
    }
    j = settings;
    settings_file << j.dump(2);
    // Parameters
    fs::path param_path = dir/"parameters.json";
    std::ofstream param_file{param_path};
    if(!param_file){
        std::cerr << "Can not open file at \"" << param_path << "\" for writing parameter sets";
        success = false;
    }
    j = json{};
    for(const auto& pair: params){
        const auto& plugin = pair.first;
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
    if(!preset_file){
        std::cerr << "Can not open file at \"" << preset_path << "\" for writing IO-presets";
        success = false;
    }
    j = json{};
    for(const auto& pair: presets){
        const auto& plugin = pair.first;
        if(!(plugin->arguments & IO::Plugin::Preset)) continue;
        const auto& com = plugin->command;
        j[com] = json{};
        auto& tmp = j[com];
        for(const auto& preset: pair.second){
            tmp[preset.first] = *preset.second;
        }
    }
    preset_file << j.dump(2);
    return;
}
