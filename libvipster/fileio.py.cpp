#include <pybind11/functional.h>
#include "nlohmann/json.hpp"
#include "pyvipster.h"
#include "fileio.h"
#include "configfile.h"
#include <map>
#include <pybind11/stl.h>

namespace Vipster::IO {
std::ostream& operator<<(std::ostream& os, const Plugin *p){
    return os << "Plugins."+p->command;
}

// TODO: nicer representation, json maybe misleading in a python context
std::ostream& operator<<(std::ostream& os, const Parameter &p){
    nlohmann::json j;
    to_json(j, p);
    return os << j;
}

std::ostream& operator<<(std::ostream& os, const Preset &p){
    nlohmann::json j;
    to_json(j, p);
    return os << j;
}

}

PYBIND11_MAKE_OPAQUE(Vipster::IO::Parameters);
PYBIND11_MAKE_OPAQUE(Vipster::IO::Parameters::mapped_type);
PYBIND11_MAKE_OPAQUE(Vipster::IO::Presets);
PYBIND11_MAKE_OPAQUE(Vipster::IO::Presets::mapped_type);

namespace Vipster::Py{

void IO(py::module& m, const ConfigState& state, bool enableRead){
    auto io = m.def_submodule("IO");

    py::class_<IO::Data>(io, "Data")
        .def(py::init<>())
        .def_readwrite("mol", &IO::Data::mol)
        .def_readwrite("param", &IO::Data::param)
//        .def_readwrite("data", &IO::Data::data)
    ;

    if(enableRead){
        /*
         * read a file
         */
        m.def("readFile",[&state](std::string fn){
            auto data = readFile(fn, std::get<2>(state));
            if(data.data.empty()){
                return py::make_tuple(data.mol, std::move(data.param), py::none());
            }else{
                py::list l{};
                for(auto& d: data.data){
                    l.append(d.release());
                }
                return py::make_tuple(data.mol, std::move(data.param), l);
            }
        }, "filename"_a);
        m.def("readFile",[](std::string fn, const IO::Plugin* plug){
            auto data = readFile(fn, plug);
            if(data.data.empty()){
                return py::make_tuple(data.mol, std::move(data.param), py::none());
            }else{
                py::list l{};
                for(auto& d: data.data){
                    l.append(d.release());
                }
                return py::make_tuple(data.mol, std::move(data.param), l);
            }
        }, "filename"_a, "format"_a);
    }

    /*
     * write a file
     *
     * falling back to default-preset/param
     */
    m.def("writeFile", [&state](const std::string &fn, const IO::Plugin* plug, const Molecule &m,
            std::optional<size_t> idx={},
            std::optional<IO::Parameter> p={},
            std::optional<IO::Preset> c={}){
        if(!p && plug->makeParam){
            p = plug->makeParam();
        }
        if(!c && plug->makePreset){
            c = plug->makePreset();
        }
        return writeFile(fn, plug, m, idx, p, c);
        },
          "filename"_a, "format"_a, "molecule"_a,
          "index"_a=std::nullopt, "param"_a=nullptr, "config"_a=nullptr);
    m.def("writeString", [&state](const IO::Plugin* plug, const Molecule &m,
            std::optional<size_t> idx={},
            std::optional<IO::Parameter> p={},
            std::optional<IO::Preset> c={}){
        if(!plug->writer){
            throw IO::Error{"Read-only format"};
        }
        if(!p && plug->makeParam){
            p = plug->makeParam();
        }
        if(!c && plug->makePreset){
            c = plug->makePreset();
        }
        if(!idx){
            idx = m.getNstep()-1;
        }
        std::stringstream ss{};
        auto success = plug->writer(m, ss, p, c, *idx);
        if(success){
            return ss.str();
        }else{
            return std::string{};
        }
    },"format"_a, "molecule"_a, "index"_a=std::nullopt, "param"_a=nullptr, "config"_a=nullptr);

    // Expose plugins
    py::class_<IO::Plugin>(io, "__Plugin")
        .def("__repr__", [](const IO::Plugin *p){return "Plugins."+p->command;})
        .def_readonly("name", &IO::Plugin::name)
        .def_readonly("extension", &IO::Plugin::extension)
        .def_readonly("command", &IO::Plugin::command)
        .def_property_readonly("hasParser", [](const IO::Plugin*const p){
            return static_cast<bool>(p->parser);})
        .def_property_readonly("hasWriter", [](const IO::Plugin*const p){
            return static_cast<bool>(p->writer);})
        .def_property_readonly("hasParameters", [](const IO::Plugin*const p){
            return static_cast<bool>(p->makeParam);})
        .def_property_readonly("hasPresets", [](const IO::Plugin*const p){
            return static_cast<bool>(p->makePreset);})
    ;

    // Expose parameters and presets
    py::bind_map<IO::Parameters::mapped_type>(io, "__StrParMap__");
    py::bind_map<IO::Parameters>(io, "__Parameters")
        .def("__repr__", [](IO::Parameters& p){
            std::ostringstream s;
            s << "__Parameters{";
            bool f = false;
            for (auto const &kv:p){
                if (f)
                    s << ", ";
                s << kv.first << ": __StrParMap__";
                f = true;
            }
            s << '}';
            return s.str();
        })
    ;

    py::bind_map<IO::Presets::mapped_type>(io, "__StrPresMap__");
    py::bind_map<IO::Presets>(io, "__IOPresets")
        .def("__repr__", [](IO::Presets& p){
            std::ostringstream s;
            s << "__Parameters{";
            bool f = false;
            for (auto const &kv:p){
                if (f)
                    s << ", ";
                s << kv.first << ": __StrPresMap__";
                f = true;
            }
            s << '}';
            return s.str();
        })
    ;

    /*
     * Initialize plugins' parameters and presets
     */
    py::class_<IO::Parameter>(io, "__BaseParam")
        .def("__repr__", [](const IO::Parameter &p){
            std::ostringstream s;
            nlohmann::json j;
            to_json(j, p);
            s << p.getFmt()->command << "-Parameters" << j;
            return s.str();
        });
    py::class_<IO::Preset>(io, "__BasePreset")
        .def("__repr__", [](const IO::Preset &p){
            std::ostringstream s;
            nlohmann::json j;
            to_json(j, p);
            s << p.getFmt()->command << "-Parameters" << j;
            return s.str();
        });

}
}
