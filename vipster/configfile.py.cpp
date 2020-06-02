#include <pybind11/eval.h>
#include <filesystem>

#include "configfile.py.h"
#include "configfile.h"
#include "io/plugin.py.h"

namespace fs = std::filesystem;

void Vipster::Py::config(py::module& m, ConfigState& state){
    m.attr("PeriodicTable") = std::get<0>(state);
    m.attr("Parameters") = py::cast(std::get<3>(state), py::return_value_policy::reference);
    m.attr("Presets") = py::cast(std::get<4>(state), py::return_value_policy::reference);
    auto& plugins = std::get<2>(state);
    auto& params = std::get<3>(state);
    auto& presets = std::get<4>(state);
    // try to parse python-based plugins
    auto pluginDir = getConfigDir()/"plugins";
    if(!fs::exists(pluginDir)) return;
    for(const auto& file: fs::directory_iterator(pluginDir)){
        if(file.path().extension() != ".py") continue;
        auto plug = Plugin::create(file);
        if(plug){
            std::cerr << "Loading Python-plugin " << file.path() << std::endl;
            plugins.push_back(plug);
            if(plug->makeParam){
                params[plug]["default"] = plug->makeParam();
            }
            if(plug->makePreset){
                presets[plug]["default"] = plug->makePreset();
            }
        }
    }
    // expose plugins
    auto plug = m.def_submodule("Plugins");
    for(const auto* p: plugins){
        plug.attr(p->command.c_str()) = p;
    }
}
