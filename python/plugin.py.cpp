#include "plugin.py.h"
#include "plugins/moltemplate.py.h"

#include <pybind11/eval.h>
#include <pybind11/stl.h>
#include <fmt/format.h>

using namespace Vipster;

std::vector<std::pair<std::string, py::object(*)()>> defaultPyPlugins(){
    return {
        {"moltemplate", &Plugins::moltemplate}
    };
}

void Vipster::Py::Plugins(py::module& m, ConfigState &state){
    // Create python-based Plugins
    for(const auto& [name, plug]: defaultPyPlugins()){
        try{
            Plugin::try_create(plug(), state);
        }catch(const py::error_already_set &e){
            std::cerr << "Can not load file format plugin " << name
                      << ":\n" <<  py::str(e.type().attr("__name__"))
                      << ": " << py::str(e.value()) << '\n' << std::endl;
        }
    }

    // register Plugins with Python
    py::class_<Vipster::Plugin>(m, "__Plugin")
    .def("__repr__", [](const Vipster::Plugin *p){return "Plugins."+p->command;})
    .def_readonly("name", &Vipster::Plugin::name)
    .def_readonly("extension", &Vipster::Plugin::extension)
    .def_readonly("command", &Vipster::Plugin::command)
    // Parser
    .def("parser", [](const Vipster::Plugin *p, const std::string &name, const std::string &content) -> py::object{
        if(!p->parser){
            throw Error{fmt::format("Plugins.{} has no parser", p->command)};
        }
        auto contentstream = std::stringstream{content};
        return py::cast(std::get<0>(p->parser(name, contentstream)));
    })
    .def_property_readonly("hasParser", [](const Vipster::Plugin*const p){
        return static_cast<bool>(p->parser);})
    // Writer
    .def("writer", [](const Vipster::Plugin *plug,
                      const Molecule& m,
                      std::optional<size_t> idx,
                      const std::optional<Parameter>& p,
                      const std::optional<Preset>& c) -> py::object{
        if(!plug->writer){
            throw Error{fmt::format("Plugins.{} has no writer", plug->command)};
        }
        if(!idx){
            idx = m.getNstep()-1;
        }
        std::stringstream file{};
        plug->writer(m, file, p, c, *idx);
        return py::cast(file.str());
    }, "molecule"_a, "index"_a=std::nullopt, "parameters"_a=std::nullopt, "preset"_a=std::nullopt)
    .def_property_readonly("hasWriter", [](const Vipster::Plugin*const p){
        return static_cast<bool>(p->writer);})
    // Parameters
    .def("makeParameters", [](const Vipster::Plugin *const p){
        if(!p->makeParam){
            throw Error{fmt::format("Plugins.{} has no parameters", p->command)};
        }
        return py::cast(p->makeParam);
    })
    .def_property_readonly("hasParameters", [](const Vipster::Plugin*const p){
        return static_cast<bool>(p->makeParam);})
    // Presets
    .def("makePreset", [](const Vipster::Plugin*const p){
        if(!p->makeParam){
            throw Error{fmt::format("Plugins.{} has no IO presets", p->command)};
        }
        return py::cast(p->makePreset());
    })
    .def_property_readonly("hasPresets", [](const Vipster::Plugin*const p){
        return static_cast<bool>(p->makePreset);})
    ;
}

std::ostream& Vipster::operator<<(std::ostream& os, const Plugin *p){
    return os << "vipster.Plugins."+p->command;
}

// implement Py::Plugin
IOTuple Py::Plugin::parser_impl(const std::string& name, std::istream &file){
    try{
        std::string filestring{std::istreambuf_iterator{file}, {}};
        auto mol = py::cast<Molecule>(pyReader(name, filestring));
        return IOTuple{mol, std::nullopt, DataList{}};
    }catch(const py::error_already_set& e){
        throw IOError{e.what()};
    }
}

bool Py::Plugin::writer_impl(const Molecule &m, std::ostream &file,
                 const std::optional<Parameter>& p,
                 const std::optional<Preset>& c,
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
        throw IOError{e.what()};
    }catch(const py::cast_error&){
        throw IOError{"Could not convert Molecule to file"};
    }
}

Parameter Py::Plugin::makeParam_impl(){
    try{
        return {this, param};
    }catch(const py::error_already_set& e){
        throw IOError{e.what()};
    }
}

Preset Py::Plugin::makePreset_impl(){
    try{
        return {this, preset};
    }catch(const py::error_already_set& e){
        throw IOError{e.what()};
    }
}

// try to create a plugin from a Python object/Module
// if a valid, complete spec is found, will register it with `state`
bool Py::Plugin::try_create(const py::object &pyplug, ConfigState &state){
    using namespace std::placeholders;
    auto plug = new Py::Plugin{};
    // Parse strings for name, extension and CLI-flag
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
    // parser and writer must be functions
    // TODO: possible to error on signature-mismatch?
    if(py::hasattr(pyplug, "parser") &&
            py::isinstance<py::function>(pyplug.attr("parser"))){
        plug->pyReader = pyplug.attr("parser");
        plug->parser = std::bind(&Py::Plugin::parser_impl, plug, _1, _2);
    }
    if(py::hasattr(pyplug, "writer") &&
            py::isinstance<py::function>(pyplug.attr("writer"))){
        plug->pyWriter = pyplug.attr("writer");
        plug->writer = std::bind(&Py::Plugin::writer_impl, plug, _1, _2, _3, _4, _5);
    }
    // Require all spec-strings and at least one function to be a useful plugin
    if(plug->name.empty() || plug->extension.empty() || plug->command.empty() ||
       (!plug->parser && !plug->writer)){
        // plugin is incomplete, abort
        delete plug;
        return false;
    }
    // register plugin with `state`'s list of plugins
    auto &plugins = std::get<2>(state);
    plugins.push_back(plug);
    // Optional: Parameter-set definition
    if(py::hasattr(pyplug, "parameters") &&
            py::isinstance<py::dict>(pyplug.attr("parameters"))){
        auto dict = py::cast<py::dict>(pyplug.attr("parameters"));
        Parameter::map_t tmp{};
        for(const auto& [key, value]: dict){
            if(py::isinstance<py::tuple>(value) &&
                    (py::cast<py::tuple>(value).size() == 2)){
                // expect (value, doc) tuple
                auto tuple = py::cast<py::tuple>(value);
                if(py::isinstance<py::str>(tuple[0])){
                    // string
                    tmp[py::str(key)] = {py::str(tuple[0]),
                                         py::str(tuple[1])};
                }else if(py::isinstance<py::list>(tuple[0])){
                    // list of strings
                    auto list = py::cast<py::list>(tuple[0]);
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
            // found a valid parameter template, register with plugin and state
            plug->param = {tmp.begin(), tmp.end()};
            plug->makeParam = std::bind(&Py::Plugin::makeParam_impl, plug);
            auto &params = std::get<3>(state);
            params[plug]["default"] = plug->makeParam();
        }
    }
    // Optional: IO-Preset definition
    if(py::hasattr(pyplug, "preset") &&
            py::isinstance<py::dict>(pyplug.attr("preset"))){
        auto dict = py::cast<py::dict>(pyplug.attr("preset"));
        Preset::map_t tmp{};
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
                    if(idx < 0 || idx >= static_cast<int>(names.size())) continue;
                    tmp[py::str(key)] = {NamedEnum{idx, names},
                                         py::str(tuple[1])};
                }
            }
        }
        if(!tmp.empty()){
            // found a valid IO-Preset template, register with plugin and state
            plug->preset = {tmp.begin(), tmp.end()};
            plug->makePreset = std::bind(&Py::Plugin::makePreset_impl, plug);
            auto &presets = std::get<4>(state);
            presets[plug]["default"] = plug->makePreset();
        }
    }
    return true;
}
