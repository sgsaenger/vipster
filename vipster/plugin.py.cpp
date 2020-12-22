#include "plugin.py.h"
#include <pybind11/eval.h>
#include <pybind11/stl.h>

using namespace Vipster;
namespace fs = std::filesystem;

// expose plugins
void Vipster::Py::Plugins(py::module& m){
    py::class_<Vipster::Plugin>(m, "__Plugin")
    .def("__repr__", [](const Vipster::Plugin *p){return "Plugins."+p->command;})
    .def_readonly("name", &Vipster::Plugin::name)
    .def_readonly("extension", &Vipster::Plugin::extension)
    .def_readonly("command", &Vipster::Plugin::command)
    .def_property_readonly("hasParser", [](const Vipster::Plugin*const p){
        return static_cast<bool>(p->parser);})
    .def_property_readonly("hasWriter", [](const Vipster::Plugin*const p){
        return static_cast<bool>(p->writer);})
    .def_property_readonly("hasParameters", [](const Vipster::Plugin*const p){
        return static_cast<bool>(p->makeParam);})
    .def_property_readonly("hasPresets", [](const Vipster::Plugin*const p){
        return static_cast<bool>(p->makePreset);})
    ;
}

std::ostream& Vipster::operator<<(std::ostream& os, const Plugin *p){
    return os << "vipster.Plugins."+p->command;
}

// implement Py::Plugin
IOTuple Py::Plugin::parser_impl(const std::string& n, std::istream &file){
    try{
        // TODO costly: copying whole file to string first
        std::string filestring;
        file >> filestring;
        IOTuple data{Molecule{n,0}, std::optional<Parameter>{}, DataList{}};
//        pyReader(py::cast(data, py::return_value_policy::reference), filestring);
        pyReader(py::cast(std::get<0>(data), py::return_value_policy::reference),
                 filestring);
        return data;
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

Py::Plugin* Py::Plugin::create(fs::path file){
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
        plug->parser = std::bind(&Py::Plugin::parser_impl, plug, _1, _2);
    }
    if(py::hasattr(pyplug, "writer") &&
            py::isinstance<py::function>(pyplug.attr("writer"))){
        plug->pyWriter = pyplug.attr("writer");
        plug->writer = std::bind(&Py::Plugin::writer_impl, plug, _1, _2, _3, _4, _5);
    }
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
            plug->param = {tmp.begin(), tmp.end()};
            plug->makeParam = std::bind(&Py::Plugin::makeParam_impl, plug);
        }
    }
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
            plug->preset = {tmp.begin(), tmp.end()};
            plug->makePreset = std::bind(&Py::Plugin::makePreset_impl, plug);
        }
    }
    if(!plug->name.empty() && !plug->extension.empty() && !plug->command.empty()){
        return plug;
    }
    // plugin is incomplete
    delete plug;
    return nullptr;
}
