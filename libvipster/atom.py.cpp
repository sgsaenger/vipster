#include "pyvipster.h"
#include "atom.h"

namespace Vipster::Py{
void Atom(py::module& m){
    auto a = py::class_<Vipster::Atom>(m, "Atom")
        .def_property("name", [](const Vipster::Atom &a)->const std::string&{return a.name;},
                      [](Vipster::Atom &a, std::string s){a.name = s;})
        .def_property("coord", [](const Vipster::Atom &a)->const Vec&{return a.coord;},
                      [](Vipster::Atom &a, Vec c){a.coord = c;})
        .def_property("properties", [](const Vipster::Atom &a)->const AtomProperties&{return a.properties;},
                      [](Vipster::Atom &a, AtomProperties bs){a.properties = bs;})
        .def(py::self == py::self)
        .def(py::self != py::self)
    ;

    py::enum_<AtomFmt>(a, "Fmt")
        .value("Bohr", AtomFmt::Bohr)
        .value("Angstrom", AtomFmt::Angstrom)
        .value("Crystal", AtomFmt::Crystal)
        .value("Alat", AtomFmt::Alat)
    ;

    py::enum_<AtomFlag>(a, "Flag")
        .value("FixX", AtomFlag::FixX)
        .value("FixY", AtomFlag::FixY)
        .value("FixZ", AtomFlag::FixZ)
        .value("Hidden", AtomFlag::Hidden)
    ;

    py::class_<AtomFlags>(a, "Flags")
        .def("__getitem__",[](const AtomFlags &bs, AtomFlag ap){
            return bs[static_cast<uint8_t>(ap)];
        })
        .def("__setitem__",[](AtomFlags &bs, AtomFlag ap, bool val){
            bs[static_cast<uint8_t>(ap)] = val;
        })
    ;

    py::class_<AtomProperties>(a, "Properties")
        .def_readwrite("charge", &AtomProperties::charge)
        .def_readwrite("forces", &AtomProperties::forces)
        .def_readwrite("flags", &AtomProperties::flags)
    ;

}
}
