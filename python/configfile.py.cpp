#include <pybind11/eval.h>
#include <filesystem>
#include <fmt/format.h>

#include "configfile.py.h"
#include "plugin.py.h"

namespace fs = std::filesystem;

void Vipster::Py::config(py::module& m, ConfigState& state){
    m.attr("PeriodicTable") = std::get<0>(state);
    m.attr("Parameters") = py::cast(std::get<3>(state), py::return_value_policy::reference);
    m.attr("Presets") = py::cast(std::get<4>(state), py::return_value_policy::reference);
    // create import function that takes absolute paths
    py::exec("import importlib.util as __util\n"
             "def __importAbsolute(name, path):\n"
             "  spec = __util.spec_from_file_location(name, path)\n"
             "  mod = __util.module_from_spec(spec)\n"
             "  spec.loader.exec_module(mod)\n"
             "  return mod"
             );
    py::object import = py::globals()["__importAbsolute"];
    // try to parse python-based plugins
    auto pluginDir = getConfigDir()/"plugins";
    if(fs::exists(pluginDir)){
        for(const fs::path& file: fs::directory_iterator(pluginDir)){
            if(file.extension() != ".py") continue;
            try{
                py::object pyplug = import(file.stem().string(), file.string());
                if(Plugin::try_create(pyplug, state)){
                    std::cerr << "Loading Python-plugin " << file << std::endl;
                }else{
                    std::cerr << "Could not load Python-plugin at " << file
                              << "\nPlugin incomplete." << std::endl;
                }
            }catch(const py::error_already_set &e){
                std::cerr << fmt::format("Could not load Python-plugin at {}:\n{}",
                                         file.string(), e.what()) << std::endl;
                continue;
            }
        }
    }
    // expose plugins
    auto plug = m.def_submodule("Plugins");
    for(const Vipster::Plugin* p: std::get<2>(state)){
        plug.attr(p->command.c_str()) = p;
    }
}
