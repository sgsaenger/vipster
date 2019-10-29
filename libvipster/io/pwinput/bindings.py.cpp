#include "pyvipster.h"
#include "plugin.h"

namespace Vipster::Py{
void PWInput(py::module& m){
    auto p = py::class_<IO::PWParam, IO::BaseParam>(m, "PWParam")
        .def_readwrite("control", &IO::PWParam::control)
        .def_readwrite("system", &IO::PWParam::system)
        .def_readwrite("electrons", &IO::PWParam::electrons)
        .def_readwrite("ions", &IO::PWParam::ions)
        .def_readwrite("cell", &IO::PWParam::cell)
    ;

    auto c = py::class_<IO::PWPreset, IO::BasePreset>(m, "PWPreset")
        .def_readwrite("atoms", &IO::PWPreset::atoms)
        .def_readwrite("cell", &IO::PWPreset::cell)
    ;
    py::enum_<IO::PWPreset::AtomFmt>(c, "AtomFmt")
        .value("Bohr", IO::PWPreset::AtomFmt::Bohr)
        .value("Angstrom", IO::PWPreset::AtomFmt::Angstrom)
        .value("Crystal", IO::PWPreset::AtomFmt::Crystal)
        .value("Alat", IO::PWPreset::AtomFmt::Alat)
//        .value("Current", IO::PWPreset::AtomFmt::Current) // TODO: makes sense to expose this?
    ;
    py::enum_<IO::PWPreset::CellFmt>(c, "CellFmt")
        .value("Angstrom", IO::PWPreset::CellFmt::Angstrom)
        .value("Bohr", IO::PWPreset::CellFmt::Bohr)
//        .value("Current", IO::PWPreset::CellFmt::Current) // TODO: makes sense to expose this?
    ;
}
}
