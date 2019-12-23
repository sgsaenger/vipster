#include <pybind11/eval.h>
#include "pyvipster.h"
#include "configfile.h"
#include "io/plugin.h"

namespace fs = std::filesystem;

namespace Vipster::Py{

class Plugin : public IO::Plugin
{
private:
    Plugin()=default;
    IO::Data parser_impl(const std::string& n, std::istream &file){
        try{
            // TODO costly: copying whole file to string first
            std::string filestring;
            file >> filestring;
            return py::cast<IO::Data>(pyReader(n, filestring));
        }catch(const py::error_already_set& e){
            throw IO::Error{e.what()};
        }catch(const py::cast_error&){
            throw IO::Error{"Could not parse a Molecule"};
        }
    }
    bool writer_impl(const Molecule &m, std::ostream &file,
                     const std::optional<IO::Parameter>& p,
                     const std::optional<IO::Preset>& c,
                     size_t idx){
        try{
            std::string filestring = py::cast<std::string>(pyWriter(m, p, c, idx));
            if(!filestring.empty()){
                file << filestring;
                return true;
            }else{
                return false;
            }
        }catch(const py::error_already_set& e){
            throw IO::Error{e.what()};
        }catch(const py::cast_error&){
            throw IO::Error{"Could not convert Molecule to file"};
        }
    }
    IO::Parameter makeParam_impl(){
        try{
            return {this, param};
        }catch(const py::error_already_set& e){
            throw IO::Error{e.what()};
        }
    }
    IO::Preset makePreset_impl(){
        try{
            return {this, preset};
        }catch(const py::error_already_set& e){
            throw IO::Error{e.what()};
        }
    }
    py::object pyReader;
    py::object pyWriter;
    IO::Parameter::BaseMap param;
    IO::Preset::BaseMap preset;
public:
    static Plugin* create(fs::path file){
        using namespace std::placeholders;
        auto plug = new Py::Plugin{};
        py::object pyplug;
        try{
            py::exec("import importlib.util as __util\n"
                     "def __importAbsolute(name, path):\n"
                     "  spec = __util.spec_from_file_location(name, path)\n"
                     "  mod = __util.module_from_spec(spec)\n"
                     "  spec.loader.exec_module(mod)\n"
                     "  return mod"
                     );
            py::object import = py::globals()["__importAbsolute"];
            pyplug = import(file.stem().string(), file.string());
        }catch(...){
            // plugin is broken
            delete plug;
            return nullptr;
        }
        if(py::hasattr(pyplug, "name") &&
                py::isinstance<py::str>(pyplug.attr("name"))){
            plug->name = py::str(pyplug.attr("name"));
        }
        if(py::hasattr(pyplug, "extension") &&
                py::isinstance<py::str>(pyplug.attr("extension"))){
            plug->extension = py::str(pyplug.attr("extension"));
        }
        if(py::hasattr(pyplug, "command") &&
                py::isinstance<py::str>(pyplug.attr("command"))){
            plug->command = py::str(pyplug.attr("command"));
        }
        if(py::hasattr(pyplug, "parser") &&
                py::isinstance<py::function>(pyplug.attr("parser"))){
            plug->pyReader = pyplug.attr("parser");
            plug->parser = std::bind(&Plugin::parser_impl, plug, _1, _2);
        }
        if(py::hasattr(pyplug, "writer") &&
                py::isinstance<py::function>(pyplug.attr("writer"))){
            plug->pyWriter = pyplug.attr("writer");
            plug->writer = std::bind(&Plugin::writer_impl, plug, _1, _2, _3, _4, _5);
        }
        if(py::hasattr(pyplug, "parameters") &&
                py::isinstance<py::dict>(pyplug.attr("parameters"))){
            auto dict = py::cast<py::dict>(pyplug.attr("parameters"));
            IO::Parameter::BaseMap::map_t tmp{};
            for(const auto& [key, value]: dict){
                if(py::isinstance<py::tuple>(value) &&
                        (py::cast<py::tuple>(value).size() == 2)){
                    // expect (value, doc) tuple
                    auto tuple = py::cast<py::tuple>(value);
                    if(py::isinstance<py::str>(tuple[0])){
                        // string
                        tmp[py::str(key)] = {py::str(tuple[0]),
                                             py::str(tuple[1])};
                    }else if(py::isinstance<py::tuple>(tuple[0])){
                        // list of strings
                        auto list = py::cast<py::tuple>(tuple[0]);
                        std::vector<std::string> vec{};
                        for(const auto& val: list){
                            vec.push_back(py::str(val));
                        }
                        tmp[py::str(key)] = {std::move(vec),
                                             py::str(tuple[1])};
                    }else if(py::isinstance<py::dict>(tuple[0])){
                        // dict of strings
                        auto dict = py::cast<py::dict>(tuple[0]);
                        std::map<std::string, std::string> map{};
                        for(const auto& [k, v]: dict){
                            map[py::str(k)] = py::str(v);
                        }
                        tmp[py::str(key)] = {std::move(map),
                                             py::str(tuple[1])};
                    }
                }
            }
            if(!tmp.empty()){
                plug->param = {tmp.begin(), tmp.end()};
                plug->makeParam = std::bind(&Plugin::makeParam_impl, plug);
            }
        }
        if(py::hasattr(pyplug, "preset") &&
                py::isinstance<py::dict>(pyplug.attr("preset"))){
            auto dict = py::cast<py::dict>(pyplug.attr("preset"));
            IO::Preset::BaseMap::map_t tmp{};
            for(const auto& [key, value]: dict){
                if(py::isinstance<py::tuple>(value) &&
                        (py::cast<py::tuple>(value).size() == 2)){
                    // expect (value, doc) tuple
                    auto tuple = py::cast<py::tuple>(value);
                    if(py::isinstance<py::bool_>(tuple[0])){
                        // boolean flag
                        tmp[py::str(key)] = {py::cast<bool>(tuple[0]),
                                             py::str(tuple[1])};
                    }else if(py::isinstance<py::tuple>(tuple[0]) &&
                             py::cast<py::tuple>(tuple[0]).size() == 2){
                        // NamedEnum -> (idx, (list-of-names))
                        auto tup2 = py::cast<py::tuple>(tuple[0]);
                        if(!py::isinstance<py::int_>(tup2[0]) ||
                           !py::isinstance<py::tuple>(tup2[1])){
                            continue;
                        }
                        auto idx = py::cast<int>(tup2[0]);
                        std::vector<std::string> names{};
                        for(const auto& n: py::cast<py::tuple>(tup2[1])){
                            names.push_back(py::str(n));
                        }
                        if(idx < 0 || idx >= names.size()) continue;
                        tmp[py::str(key)] = {NamedEnum{idx, names},
                                             py::str(tuple[1])};
                    }
                }
            }
            if(!tmp.empty()){
                plug->preset = {tmp.begin(), tmp.end()};
                plug->makePreset = std::bind(&Plugin::makePreset_impl, plug);
            }
        }
        if(!plug->name.empty() && !plug->extension.empty() && !plug->command.empty()){
            return plug;
        }
        // plugin is incomplete
        delete plug;
        return nullptr;
    }
};

void config(py::module& m, ConfigState& state){
    m.attr("PeriodicTable") = std::get<0>(state);
    m.attr("Parameters") = py::cast(std::get<3>(state), py::return_value_policy::reference);
    m.attr("IOPresets") = py::cast(std::get<4>(state), py::return_value_policy::reference);
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
}
