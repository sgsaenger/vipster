#include "pyvipster.h"
#include "io/cpmdinput/plugin.h"

namespace Vipster::Py{
void CPInput(py::module& m){
    auto p = py::class_<IO::CPParam>(m, "CPParam")
        .def_readwrite("info", &IO::CPParam::info)
        .def_readwrite("cpmd", &IO::CPParam::cpmd)
        .def_readwrite("system", &IO::CPParam::system)
        .def_readwrite("pimd", &IO::CPParam::pimd)
        .def_readwrite("path", &IO::CPParam::path)
        .def_readwrite("ptddft", &IO::CPParam::ptddft)
        .def_readwrite("atoms", &IO::CPParam::atoms)
        .def_readwrite("dft", &IO::CPParam::dft)
        .def_readwrite("prop", &IO::CPParam::prop)
        .def_readwrite("resp", &IO::CPParam::resp)
        .def_readwrite("linres", &IO::CPParam::linres)
        .def_readwrite("tddft", &IO::CPParam::tddft)
        .def_readwrite("hardness", &IO::CPParam::hardness)
        .def_readwrite("classic", &IO::CPParam::classic)
        .def_readwrite("exte", &IO::CPParam::exte)
        .def_readwrite("vdw", &IO::CPParam::vdw)
        .def_readwrite("qmmm", &IO::CPParam::qmmm)
    ;

    auto c = py::class_<IO::CPConfig>(m, "CPConfig")
        .def_readwrite("fmt", &IO::CPConfig::fmt)
    ;

    py::enum_<IO::CPConfig::AtomFmt>(c, "AtomFmt")
        .value("Bohr", IO::CPConfig::AtomFmt::Bohr)
        .value("Angstrom", IO::CPConfig::AtomFmt::Angstrom)
        .value("Crystal", IO::CPConfig::AtomFmt::Crystal)
        .value("Alat", IO::CPConfig::AtomFmt::Alat)
//        .value("Current", IO::CPConfig::AtomFmt::Current) // TODO: makes sense to expose this?
    ;
}
}
