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
            return py::cast<IO::Parameter&>(pyMkParam());
        }catch(const py::error_already_set& e){
            throw IO::Error{e.what()};
        }
    }
    IO::Preset makePreset_impl(){
        try{
            return py::cast<IO::Preset&>(pyMkPreset());
        }catch(const py::error_already_set& e){
            throw IO::Error{e.what()};
        }
    }
    py::object pyReader;
    py::object pyWriter;
    py::object pyMkParam;
    py::object pyMkPreset;
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
        if(py::hasattr(pyplug, "name")){
            plug->name = py::cast<std::string>(pyplug.attr("name"));
        }
        if(py::hasattr(pyplug, "extension")){
            plug->extension = py::cast<std::string>(pyplug.attr("extension"));
        }
        if(py::hasattr(pyplug, "command")){
            plug->command = py::cast<std::string>(pyplug.attr("command"));
        }
        if(py::hasattr(pyplug, "parser")){
            plug->pyReader = pyplug.attr("parser");
            plug->parser = std::bind(&Plugin::parser_impl, plug, _1, _2);
        }
        if(py::hasattr(pyplug, "writer")){
            plug->pyWriter = pyplug.attr("writer");
            plug->writer = std::bind(&Plugin::writer_impl, plug, _1, _2, _3, _4, _5);
        }
        // TODO: enable when sufficient implementation is available
//        if(py::hasattr(pyplug, "makeParam")){
//            plug->pyMkParam = pyplug.attr("makeParam");
//            plug->makeParam = std::bind(&Plugin::makeParam_impl, plug);
//        }
//        if(py::hasattr(pyplug, "makePreset")){
//            plug->pyMkPreset = pyplug.attr("makePreset");
//            plug->makePreset = std::bind(&Plugin::makePreset_impl, plug);
//        }
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
    // try to parse python-based plugins
    auto pluginDir = getConfigDir()/"plugins";
    if(!fs::exists(pluginDir)) return;
    for(const auto& file: fs::directory_iterator(pluginDir)){
        if(file.path().extension() != ".py") continue;
        auto plug = Plugin::create(file);
        if(plug){
            std::cerr << "Loading Python-plugin " << file.path() << std::endl;
            plugins.push_back(plug);
        }
    }
    // expose plugins
    auto plug = m.def_submodule("Plugins");
    for(const auto* p: plugins){
        plug.attr(p->command.c_str()) = p;
    }
}
}
