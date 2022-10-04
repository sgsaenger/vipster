#include "presets.py.h"
#include "plugin.py.h"
#include "vipster/presets.h"
#include <fmt/format.h>

using namespace Vipster;

namespace std{
std::ostream& operator<<(std::ostream &os, const Preset::mapped_type &p){
    switch(p.first.index()){
    case Preset::i_bool:
        os << (std::get<Preset::i_bool>(p.first) ? "True" : "False");
        break;
    case Preset::i_enum:
        os << '"' << std::get<Preset::i_enum>(p.first).name() << '"';
        break;
    }
    return os;
}
}

namespace Vipster{
std::ostream& operator<<(std::ostream& os, const Preset &p){
    os << p.getFmt()->command << "-Preset{";
    bool f{false};
    for(auto [key, val]: p){
        if (f) os << ", ";
        os << '"' << key << "\": " << val;
        f = true;
    }
    os << '}';
    return os;
}
std::ostream& operator<<(std::ostream &os, const PresetMap::mapped_type &p){
    os << "__StrPresMap{";
    bool f{false};
    for(auto [key, val]: p){
        if (f) os << ", ";
        os << '"' << key << "\": " << val;
        f = true;
    }
    os << '}';
    return os;
}
}


void Vipster::Py::Presets(py::module &m){
    py::class_<Preset>(m, "__Preset")
        .def("__repr__", [](const Preset &p){
                 std::ostringstream s;
                 s << p;
                 return s.str();
             })
        .def("__iter__",
             [](Preset &p){ return py::make_key_iterator(p.begin(), p.end()); },
             py::keep_alive<0,1>())
        .def("items",
             [](Preset &p){ return py::make_iterator(p.begin(), p.end()); },
             py::keep_alive<0,1>())
        .def("__getitem__",
             [](Preset &p, const std::string &k) {
                 auto it = p.find(k);
                 if (it == p.end())
                     throw py::key_error();
                 switch(it->second.first.index()){
                 case Preset::i_bool:
                     return py::cast(std::get<Preset::i_bool>(it->second.first), py::return_value_policy::reference);
                 case Preset::i_enum:
                     return py::cast(std::get<Preset::i_enum>(it->second.first).name(), py::return_value_policy::reference);
                 default:
                     throw py::value_error("Invalid enum value");
                 }
             }, py::return_value_policy::reference_internal)
        .def("__setitem__",
             [](Preset &p, const std::string &k, py::object &val){
                 auto it = p.find(k);
                 if (it == p.end())
                     throw py::key_error();
                 auto& tgt = it->second.first;
                 switch(tgt.index()){
                 case Preset::i_bool:
                     if(!py::isinstance<py::bool_>(val))
                         throw py::type_error(std::string(py::str(val.get_type()))+
                                              " can not be interpreted as a bool");
                     std::get<Preset::i_bool>(tgt) = py::bool_(val);
                     break;
                 case Preset::i_enum:
                 {
                     auto &enum_ = std::get<Preset::i_enum>(tgt);
                     if(!py::isinstance<py::str>(val) ||
                        std::find_if(enum_.begin(), enum_.end(), [val](const auto &e){return e.second==std::string(py::str(val));})
                             == enum_.end())
                     {
                         throw py::type_error(fmt::format("{} is not one of the valid values {{'{}'}}",
                                                          std::string(py::str(val)),
                                                          fmt::join(enum_.names(),"', '")));
                     }
                     std::get<Preset::i_enum>(tgt) = py::str(val);
                     break;
                 }
                 }
             })
        .def("__contains__",
             [](Preset &p, const std::string &k){
                 auto it = p.find(k);
                 return it != p.end();
             })
        .def("__len__", &Preset::size)
        .def("doc", [](Preset &p, const std::string &k){
                 auto it = p.find(k);
                 if (it == p.end())
                     throw py::key_error();
                 switch(it->second.first.index()){
                 case Preset::i_bool:
                     return fmt::format("{}:\nbool\n\n{}", it->first, it->second.second);
                 case Preset::i_enum:
                     return fmt::format("{}:\n{{'{}'}}\n\n{}",
                                        it->first,
                                        fmt::join(std::get<Preset::i_enum>(it->second.first).names(), "', '"),
                                        it->second.second);
                 default:
                     throw py::value_error("Invalid enum value");
                 }
             })
    ;

    py::bind_map<PresetMap>(m, "__PresetMap");
    py::bind_map<PresetMap::mapped_type>(m, "__StrPresMap");
}
