#include "pyvipster.h"
#include "fileio.h"
#include "configfile.h"
#include <map>

namespace Vipster::IO {
std::ostream& operator<<(std::ostream& os, const Plugin *p){
    return os << "Plugins."+p->command;
}

// TODO: nicer representation, json maybe misleading in a python context
std::ostream& operator<<(std::ostream& os, const std::unique_ptr<BaseParam> &p){
    return os << p->toJson();
}

std::ostream& operator<<(std::ostream& os, const std::unique_ptr<BasePreset> &p){
    return os << p->toJson();
}

}

/*
 * Custom map-binding routine for map containing unique_ptr
 *
 * Doesn't fit assumptions made by pybind11,
 * probably because unique_ptr is default holder_type and clashes
 *
 * Borrowed and modified from bind_map @ pybind11/stl_bind.h
 * and make_iterator @ pybind11.h
 */
template <typename Map>
auto bind_map_own(py::handle scope, const std::string &name) {
    using KeyType = typename Map::key_type;
    using holder_type = std::shared_ptr<Map>;
    using Class_ = py::class_<Map, holder_type>;


    Class_ cl(scope, name.c_str());

    cl.def("__repr__",
        [name](const Map &m) {
            std::ostringstream s;
            s << name << '{';
            bool f = false;
            for (auto const &kv : m) {
                if (f) s << ", ";
                s << '"' << kv.first << "\": " << kv.second;
                f = true;
            }
            s << '}';
            return s.str();
        }
    );

    cl.def("__bool__",
        [](const Map &m) -> bool { return !m.empty(); }
    );

    cl.def("__iter__",
           [](Map &m) { return py::make_key_iterator(m.begin(), m.end()); },
           py::keep_alive<0, 1>() /* Essential: keep list alive while iterator exists */
    );

    cl.def("items",
           [](Map &m) {
               using Iterator = decltype(std::declval<Map>().begin());
               using Sentinel = decltype(std::declval<Map>().end());
               constexpr auto Policy = py::return_value_policy::reference_internal;
               using state = py::detail::iterator_state<Iterator, Sentinel,
                                false, Policy>;
               if(!py::detail::get_type_info(typeid(state), false)){
                   py::class_<state>(py::handle(), "iterator", pybind11::module_local())
                       .def("__iter__", [](state &s) -> state& { return s; })
                       .def("__next__", [](state &s) -> typename Map::iterator& {
                           if (!s.first_or_done)
                               ++s.it;
                           else
                               s.first_or_done = false;
                           if (s.it == s.end) {
                               s.first_or_done = true;
                               throw py::stop_iteration();
                           }
                           return s.it;
                       }, Policy);
               }
               return py::cast(state{m.begin(), m.end(), true});
           },
           py::keep_alive<0, 1>() /* Essential: keep list alive while iterator exists */
    );

    cl.def("__getitem__",
        [](Map &m, const KeyType &k){
            auto it = m.find(k);
            if (it == m.end())
              throw py::key_error();
           return it->second.get();
        },
        py::return_value_policy::reference_internal
    );

    cl.def("__contains__",
        [](Map &m, const KeyType &k) -> bool {
            auto it = m.find(k);
            if (it == m.end())
              return false;
           return true;
        }
    );

    cl.def("__len__", &Map::size);

    return cl;
}


namespace Vipster::Py{
void PWInput(py::module&);
void LmpInput(py::module&);
void XYZ(py::module&);
void CPInput(py::module&);
void ORCA(py::module&);
void POSCAR(py::module&);

void IO(py::module& m, const ConfigState& state, bool enableRead){
    auto io = m.def_submodule("IO");

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

    /*
     * write a file
     *
     * falling back to default-preset/param
     */
    m.def("writeFile", [&state](const std::string &fn, const IO::Plugin* plug, const Molecule &m,
                       size_t idx=-1ul, const IO::BaseParam *p=nullptr, const IO::BasePreset *c=nullptr){
        if(!p && plug->makeParam){
            p = std::get<3>(state).at(plug).at("default").get();
        }
        if(!c && plug->makePreset){
            c = std::get<4>(state).at(plug).at("default").get();
        }
        return writeFile(fn, plug, m, idx, p, c);
        },
          "filename"_a, "format"_a, "molecule"_a,
          "index"_a=-1ul, "param"_a=nullptr, "config"_a=nullptr);

    // Expose plugins
    py::class_<IO::Plugin>(io, "__Plugin")
        .def("__repr__", [](const IO::Plugin *p){return "Plugins."+p->command;})
        .def_readonly("name", &IO::Plugin::name)
        .def_readonly("extension", &IO::Plugin::extension)
        .def_readonly("command", &IO::Plugin::command)
        //TODO: shall be removed on C++ side, how can it be replaced here?
        .def_readonly("arguments", &IO::Plugin::arguments)
    ;

    /*
     * Expose parameters and presets
     *
     * Circumvent bind_map for inner map
     */
    bind_map_own<IO::Parameters::mapped_type>(io, "__StrParMap__");
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

    bind_map_own<IO::Presets::mapped_type>(io, "__StrPresMap__");
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
    py::class_<IO::BaseParam>(io, "__BaseParam")
        .def("__repr__", [](const IO::BaseParam &p){
            std::ostringstream s;
            s << p.getFmt()->command << "-Parameters" << p.toJson();
            return s.str();
        });
    py::class_<IO::BasePreset>(io, "__BasePreset")
        .def("__repr__", [](const IO::BasePreset &p){
            std::ostringstream s;
            s << p.getFmt()->command << "-IOPreset" << p.toJson();
            return s.str();
        });
    PWInput(io);
    LmpInput(io);
    XYZ(io);
    CPInput(io);
    ORCA(io);
    POSCAR(io);

}
}
