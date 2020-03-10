#include "pyvipster.h"
#include "atom.h"

namespace Vipster::Py{
void Atom(py::module& m){
    py::enum_<AtomFmt>(m, "Fmt")
        .value("Crystal", AtomFmt::Crystal)
        .value("Alat", AtomFmt::Alat)
        .value("Angstrom", AtomFmt::Angstrom)
        .value("Bohr", AtomFmt::Bohr)
    ;

    py::enum_<AtomFlag>(m, "Flag")
        .value("FixX", AtomFlag::FixX)
        .value("FixY", AtomFlag::FixY)
        .value("FixZ", AtomFlag::FixZ)
        .value("Hidden", AtomFlag::Hidden)
    ;

    py::class_<AtomFlags>(m, "Flags")
        .def("__getitem__",[](const AtomFlags &bs, AtomFlag ap){
            return static_cast<bool>(bs[static_cast<uint8_t>(ap)]);
        })
        .def("__setitem__",[](AtomFlags &bs, AtomFlag ap, bool val){
            bs[static_cast<uint8_t>(ap)] = val;
        })
    ;

    py::class_<AtomProperties>(m, "Properties")
        .def_readwrite("charge", &AtomProperties::charge)
        .def_readwrite("forces", &AtomProperties::forces)
        .def_readwrite("flags", &AtomProperties::flags)
    ;
}
}
