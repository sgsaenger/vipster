#include "pyvipster.h"
#include "fileio.h"
#include <map>

namespace std{
ostream& operator<<(ostream& os, const Vipster::IO::Plugin* p){
    return os << p->command;
}
//TODO: what's wrong with map_if_insertion_operator?
//ostream& operator<<(ostream& os, const Vipster::IO::BaseParam &p){
//    return os << p.toJson();
//}
//ostream& operator<<(ostream& os, const Vipster::IO::BaseParam *p){
//    return os << p->toJson();
//}
//ostream& operator<<(ostream& os, const std::unique_ptr<Vipster::IO::BaseParam> &p){
//    return os << p->toJson();
//}
}

namespace Vipster::Py{
void PWInput(py::module&);
void LmpInput(py::module&);
void XYZ(py::module&);
void CPInput(py::module&);
void ORCA(py::module&);
void POSCAR(py::module&);

void IO(py::module& m){
    auto io = m.def_submodule("IO");

    // TODO: state handling
    m.def("readFile",[](std::string fn){
        auto data = readFile(fn);
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
    m.def("readFile",[](std::string fn, std::string fmt){
        const auto plugins = IO::defaultPlugins();
        auto pos = std::find_if(plugins.begin(), plugins.end(),
            [&](const auto& plug){return plug->command == fmt;});
        if(pos == plugins.end()){
            throw IO::Error{"Invalid format given: "+fmt};
        }
        auto data = readFile(fn, *pos);
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
     * TODO: provide wrapper
     *
     * fall back to default-preset/param
     */
    m.def("writeFile", &writeFile, "filename"_a, "format"_a, "molecule"_a,
          "param"_a=nullptr, "config"_a=nullptr, "index"_a=-1ul);

    // Expose plugins
    py::class_<IO::Plugin>(io, "__Plugin")
            // TODO: this fails on temporary objects!
            // works when wrapper is bound to a name
        .def(py::init([](const std::string &s){
            const auto plugins = IO::defaultPlugins();
            auto pos = std::find_if(plugins.begin(), plugins.end(),
                [&](const auto &plug){return plug->command == s;});
            if(pos == plugins.end()){
                throw IO::Error{"Invalid format given: "+s};
            }
            return const_cast<IO::Plugin*>(*pos);
        }), py::return_value_policy::reference_internal)
        .def("__repr__", [](const IO::Plugin *p){return p->command;})
    ;
    py::bind_vector<IO::Plugins>(io, "__Plugins");

    /*
     * Expose parameters
     *
     * Circumvent bind_map for inner map
     */
    using _TPM = IO::Parameters::mapped_type;
    auto cl = py::class_<_TPM>(m, "__StrParMap__")
        .def("__getitem__", [](const _TPM &m, const std::string &s){
            auto tmp = m.find(s);
            if(tmp != m.end()){
                return tmp->second->copy();
            }else{
                throw py::key_error();
            }})
        .def("__setitem__", [](_TPM &m, const std::string &s,
                               const IO::BaseParam *p){
                m.insert_or_assign(s, p->copy());
            })
        .def("__repr__",
            [](_TPM &m) {
            std::ostringstream s;
            s << "__StrParMap__{";
            bool f = false;
            for (auto const &kv : m) {
                if (f)
                    s << ", ";
                s << kv.first << ": " << kv.second->toJson();
                f = true;
            }
            s << '}';
            return s.str();
        }, "Return the canonical string representation of this map.")
    ;
    py::bind_map<IO::Parameters>(io, "__Parameters")
//        .def("__getitem__", [](const IO::Parameters ){})
    ;

//    {
//        auto cl = py::class_<IO::Parameters::mapped_type>(io, "__StrParMap");
//            cl.def(py::init<>());
//        using Map = IO::Parameters::mapped_type;
//        using KeyType = Map::key_type;
//        using MappedType = Map::mapped_type;
//        using Class_ = py::class_<Map, std::unique_ptr<Map>>;
//        using namespace py;
//        // Register stream insertion operator (if possible)
//        detail::map_if_insertion_operator<Map, Class_>(cl, "__StrParMap");
//        cl.def("__bool__",
//            [](const Map &m) -> bool { return !m.empty(); },
//            "Check whether the map is nonempty"
//        );
//        cl.def("__iter__",
//               [](Map &m) { return make_key_iterator(m.begin(), m.end()); },
//               keep_alive<0, 1>() /* Essential: keep list alive while iterator exists */
//        );
////        cl.def("items",
////               [](Map &m) { return make_iterator(m.begin(), m.end()); },
////               keep_alive<0, 1>() /* Essential: keep list alive while iterator exists */
////        );
//        cl.def("__getitem__",
//            [](Map &m, const KeyType &k) -> MappedType & {
//                auto it = m.find(k);
//                if (it == m.end())
//                  throw key_error();
//               return it->second;
//            },
//            return_value_policy::reference_internal // ref + keepalive
//        );
//        cl.def("__contains__",
//            [](Map &m, const KeyType &k) -> bool {
//                auto it = m.find(k);
//                if (it == m.end())
//                  return false;
//               return true;
//            }
//        );
//        // Assignment provided only if the type is copyable
//        detail::map_assignment<Map, Class_>(cl);
//        cl.def("__delitem__",
//               [](Map &m, const KeyType &k) {
//                   auto it = m.find(k);
//                   if (it == m.end())
//                       throw key_error();
//                   m.erase(it);
//               }
//        );
//        cl.def("__len__", &Map::size);
//    }



    py::bind_map<IO::Presets>(io, "__IOPresets");

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

    /*
     * Initialize plugins
     */
    PWInput(io);
    LmpInput(io);
    XYZ(io);
    CPInput(io);
    ORCA(io);
    POSCAR(io);

}
}
