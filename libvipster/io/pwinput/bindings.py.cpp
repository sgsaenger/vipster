#include "pyvipster.h"
#include "io/pwinput/plugin.h"

namespace Vipster::Py{
void PWInput(py::module& m){
    auto p = py::class_<IO::PWParam, IO::BaseParam>(m, "PWParam")
        .def_readwrite("control", &IO::PWParam::control)
        .def_readwrite("system", &IO::PWParam::system)
        .def_readwrite("electrons", &IO::PWParam::electrons)
        .def_readwrite("ions", &IO::PWParam::ions)
        .def_readwrite("cell", &IO::PWParam::cell)
    ;

    auto c = py::class_<IO::PWConfig, IO::BaseConfig>(m, "PWConfig")
        .def_readwrite("atoms", &IO::PWConfig::atoms)
        .def_readwrite("cell", &IO::PWConfig::cell)
    ;
    py::enum_<IO::PWConfig::AtomFmt>(c, "AtomFmt")
        .value("Bohr", IO::PWConfig::AtomFmt::Bohr)
        .value("Angstrom", IO::PWConfig::AtomFmt::Angstrom)
        .value("Crystal", IO::PWConfig::AtomFmt::Crystal)
        .value("Alat", IO::PWConfig::AtomFmt::Alat)
//        .value("Current", IO::PWConfig::AtomFmt::Current) // TODO: makes sense to expose this?
    ;
    py::enum_<IO::PWConfig::CellFmt>(c, "CellFmt")
        .value("Angstrom", IO::PWConfig::CellFmt::Angstrom)
        .value("Bohr", IO::PWConfig::CellFmt::Bohr)
//        .value("Current", IO::PWConfig::CellFmt::Current) // TODO: makes sense to expose this?
    ;
}
}
