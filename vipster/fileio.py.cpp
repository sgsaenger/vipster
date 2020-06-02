#include <map>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <nlohmann/json.hpp>
#include "fileio.py.h"
#include "fileio.h"
#include "io/plugin.py.h"
#include "io/presets.py.h"
#include "io/parameters.py.h"

void Vipster::Py::IO(py::module& m, const ConfigState& state, bool enableRead){
    auto io = m.def_submodule("IO");

    py::class_<IO::Data>(io, "Data")
        .def(py::init<>())
        .def_readwrite("mol", &IO::Data::mol)
        .def_readwrite("param", &IO::Data::param)
//        .def_readwrite("data", &IO::Data::data) // TODO
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
    m.def("writeFile", [](const std::string &fn, const IO::Plugin* plug, const Molecule &m,
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
    m.def("writeString", [](const IO::Plugin* plug, const Molecule &m,
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

    // Expose parameters and presets
    Py::Plugins(io);
    Py::Parameters(io);
    Py::Presets(io);
}
