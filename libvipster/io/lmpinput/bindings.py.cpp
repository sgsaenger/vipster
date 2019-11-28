#include "pyvipster.h"
#include "plugin.h"

namespace Vipster::Py{
void LmpInput(py::module &m){
    auto c = py::class_<IO::LmpPreset, IO::BasePreset>(m, "LmpPreset");

    py::enum_<IO::LmpPreset::AtomStyle>(c, "AtomStyle")
        .value("Angle", IO::LmpPreset::AtomStyle::Angle)
        .value("Atomic", IO::LmpPreset::AtomStyle::Atomic)
        .value("Bond", IO::LmpPreset::AtomStyle::Bond)
        .value("Charge", IO::LmpPreset::AtomStyle::Charge)
        .value("Full", IO::LmpPreset::AtomStyle::Full)
        .value("Molecular", IO::LmpPreset::AtomStyle::Molecular)
    ;

    c.def(py::init<IO::LmpPreset::AtomStyle, bool, bool, bool, bool>(),
          "style"_a=IO::LmpPreset::AtomStyle::Atomic,
          "bonds"_a=false, "angles"_a=false, "dihedrals"_a=false, "impropers"_a=false)
      .def_readwrite("style", &IO::LmpPreset::style)
      .def_readwrite("angles", &IO::LmpPreset::angles)
      .def_readwrite("bonds", &IO::LmpPreset::bonds)
      .def_readwrite("dihedrals", &IO::LmpPreset::dihedrals)
      .def_readwrite("impropers", &IO::LmpPreset::impropers)
    ;
}
}
