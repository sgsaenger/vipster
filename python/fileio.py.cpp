#include <map>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <fmt/format.h>
#include "fileio.py.h"
#include "vipster/fileio.h"

void Vipster::Py::FileIO(py::module& m, const ConfigState& state, bool enableRead){
    if(enableRead){
        /*
         * read a file
         */
        m.def("readFile",[&state](std::string fn){
            if(const auto plug = guessFmt(fn, std::get<2>(state))){
                return readFile(fn, plug);
            }else{
                throw IOError{fmt::format("Could not deduce format of file \"{}\""
                                          "\nPlease specify format explicitely", fn)};
            }
        }, "filename"_a);
        m.def("readFile",[](std::string fn, const Plugin* plug){
            return readFile(fn, plug);
        }, "filename"_a, "format"_a);
    }

    /*
     * write a file
     *
     * falling back to default-preset/param
     */
    m.def("writeFile", [](const std::string &fn, const Plugin* plug, const Molecule &m,
            std::optional<size_t> idx={},
            std::optional<Parameter> p={},
            std::optional<Preset> c={}){
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
    m.def("writeString", [](const Plugin* plug, const Molecule &m,
            std::optional<size_t> idx={},
            std::optional<Parameter> p={},
            std::optional<Preset> c={}){
        if(!plug->writer){
            throw IOError{"Read-only format"};
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
}
