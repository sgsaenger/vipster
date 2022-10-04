#include "parameters.py.h"
#include "plugin.py.h"
#include "vipster/parameters.h"
#include <fmt/format.h>

using namespace Vipster;

namespace std{
std::ostream& operator<<(std::ostream &os, const Parameter::mapped_type &p){
    bool f{false};
    switch (p.first.index()){
    case Parameter::i_str:
        os << '"' << std::get<Parameter::i_str>(p.first) << '"';
        break;
    case Parameter::i_strvec:
        os << '[';
        for(const auto &v: std::get<Parameter::i_strvec>(p.first)){
            if(f) os << ", ";
            os << '"' << v << '"';
            f = true;
        }
        os << ']';
        break;
    case Parameter::i_strmap:
        os << '{';
        for(auto [key, val]: std::get<Parameter::i_strmap>(p.first)){
            if(f) os << ", ";
            os << '"' << key << "\": \"" << val << '"';
            f = true;
        }
        os << '}';
        break;
    }
    return os;
}
}
namespace Vipster {
std::ostream& operator<<(std::ostream &os, const Parameter &p){
    os << p.getFmt()->command << "-Parameter{";
    bool f{false};
    for(auto [key, val]: p){
        if (f) os << ", ";
        os << '"' << key << "\": " << val;
        f = true;
    }
    os << '}';
    return os;
}
std::ostream& operator<<(std::ostream &os, const ParameterMap::mapped_type &p){
    os << "__StrParMap{";
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

void Py::Parameters(py::module &m){
    py::class_<Parameter>(m, "__Parameter")
        .def("__repr__", [](const Parameter &p){
            std::ostringstream s;
            s << p;
            return s.str();
        })
        .def("__iter__",
             [](Parameter &p){ return py::make_key_iterator(p.begin(), p.end()); },
             py::keep_alive<0,1>())
        .def("items",
             [](Parameter &p){ return py::make_iterator(p.begin(), p.end()); },
             py::keep_alive<0,1>())
        .def("__getitem__",
             [](Parameter &p, const std::string &k) {
                 auto it = p.find(k);
                 if (it == p.end())
                     throw py::key_error();
                 switch(it->second.first.index()){
                 case Parameter::i_str:
                     return py::cast(std::get<Parameter::i_str>(it->second.first), py::return_value_policy::reference);
                 case Parameter::i_strvec:
                     return py::cast(std::get<Parameter::i_strvec>(it->second.first), py::return_value_policy::reference);
                 case Parameter::i_strmap:
                     return py::cast(std::get<Parameter::i_strmap>(it->second.first), py::return_value_policy::reference);
                 default:
                     throw py::value_error("Invalid enum value");
                 }
             }, py::return_value_policy::reference_internal)
        .def("__setitem__",
             [](Parameter &p, const std::string &k, py::object &val){
                 auto it = p.find(k);
                 if (it == p.end())
                     throw py::key_error();
                 auto& tgt = it->second.first;
                 size_t idx{0};
                 switch(tgt.index()){
                 case Parameter::i_str:
                     if(!py::isinstance<py::str>(val))
                         throw py::value_error(std::string(py::str(val.get_type()))+
                                               " can not be interpreted as a string");
                     std::get<Parameter::i_str>(tgt) = py::str(val);
                     break;
                 case Parameter::i_strvec:
                     if(!py::isinstance<py::iterable>(val))
                         throw py::value_error(std::string(py::str(val.get_type()))+
                                               " can not be interpreted as a list");
                     for(const auto &v: py::iterable(val)){
                         std::get<Parameter::i_strvec>(tgt)[idx++] = py::str(v);
                     }
                     break;
                 case Parameter::i_strmap:
                     if(!py::isinstance<py::dict>(val))
                         throw py::value_error(std::string(py::str(val.get_type()))+
                                               " can not be interpreted as a dict");
                     for(auto [key, v]: py::dict(val)){
                         std::get<Parameter::i_strmap>(tgt)[py::str(key)] = py::str(v);
                     }
                     break;
                 }
             })
        .def("__contains__",
             [](Parameter &p, const std::string &k){
                 auto it = p.find(k);
                 return it != p.end();
             })
        .def("__len__", &Parameter::size)
        .def("doc", [](Parameter &p, const std::string &k){
                 auto it = p.find(k);
                 if (it == p.end())
                     throw py::key_error();
                 constexpr const char* typenames[] = {
                     "str", "[str]", "{str: str}"
                 };
                 return fmt::format("{}:\n{}\n\n{}", it->first, typenames[it->second.first.index()], it->second.second);
             })
    ;

    py::bind_map<ParameterMap>(m, "__ParameterMap");
    py::bind_map<ParameterMap::mapped_type>(m, "__StrParMap");
}
