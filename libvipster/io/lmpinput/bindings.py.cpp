#include "pyvipster.h"
#include "io/lmpinput/plugin.h"

namespace Vipster::Py{
void LmpInput(py::module &m){
    auto c = py::class_<IO::LmpConfig, IO::BaseConfig>(m, "LmpConfig")
        .def_readwrite("style", &IO::LmpConfig::style)
        .def_readwrite("angles", &IO::LmpConfig::angles)
        .def_readwrite("bonds", &IO::LmpConfig::bonds)
        .def_readwrite("dihedrals", &IO::LmpConfig::dihedrals)
        .def_readwrite("impropers", &IO::LmpConfig::impropers)
    ;

    py::enum_<IO::LmpConfig::AtomStyle>(c, "AtomStyle")
        .value("Angle", IO::LmpConfig::AtomStyle::Angle)
        .value("Atomic", IO::LmpConfig::AtomStyle::Atomic)
        .value("Bond", IO::LmpConfig::AtomStyle::Bond)
        .value("Charge", IO::LmpConfig::AtomStyle::Charge)
        .value("Full", IO::LmpConfig::AtomStyle::Full)
        .value("Molecular", IO::LmpConfig::AtomStyle::Molecular)
    ;
}
}
