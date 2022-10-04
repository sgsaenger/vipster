#include "atom.py.h"
#include "vipster/atom.h"

void Vipster::Py::Atom(py::module& m){
    py::enum_<AtomFmt>(m, "Fmt")
        .value("Crystal", AtomFmt::Crystal)
        .value("Alat", AtomFmt::Alat)
        .value("Angstrom", AtomFmt::Angstrom)
        .value("Bohr", AtomFmt::Bohr)
    ;

    auto p = py::class_<AtomProperties>(m, "Properties")
        .def_readwrite("charge", &AtomProperties::charge)
        .def_readwrite("forces", &AtomProperties::forces)
        .def_readwrite("flags", &AtomProperties::flags)
    ;

    py::enum_<AtomProperties::Flag>(p, "Flag")
        .value("FixX", AtomProperties::FixX)
        .value("FixY", AtomProperties::FixY)
        .value("FixZ", AtomProperties::FixZ)
        .value("Hidden", AtomProperties::Hidden)
    ;

    py::class_<AtomProperties::Flags>(p, "__Flags")
        .def("__getitem__",[](const AtomProperties::Flags &bs, AtomProperties::Flag ap){
            return static_cast<bool>(bs[static_cast<uint8_t>(ap)]);
        })
        .def("__setitem__",[](AtomProperties::Flags &bs, AtomProperties::Flag ap, bool val){
            bs[static_cast<uint8_t>(ap)] = val;
        })
    ;
}
