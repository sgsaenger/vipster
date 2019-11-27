#include "pyvipster.h"
#include "plugin.h"

namespace Vipster::Py{
void PWInput(py::module& m){
    auto p = py::class_<IO::PWParam, IO::BaseParam>(m, "PWParam")
        .def(py::init<IO::PWParam::Namelist,
                      IO::PWParam::Namelist,
                      IO::PWParam::Namelist,
                      IO::PWParam::Namelist,
                      IO::PWParam::Namelist,
                      std::string, std::string>(),
             "control"_a=IO::PWParam::Namelist{},
             "system"_a=IO::PWParam::Namelist{},
             "electrons"_a=IO::PWParam::Namelist{},
             "ions"_a=IO::PWParam::Namelist{},
             "cell"_a=IO::PWParam::Namelist{},
             "PPPrefix"_a=std::string{}, "PPSuffix"_a=std::string{})
        .def_readwrite("control", &IO::PWParam::control)
        .def_readwrite("system", &IO::PWParam::system)
        .def_readwrite("electrons", &IO::PWParam::electrons)
        .def_readwrite("ions", &IO::PWParam::ions)
        .def_readwrite("cell", &IO::PWParam::cell)
        .def_readwrite("PPPrefix", &IO::PWParam::PPPrefix)
        .def_readwrite("PPSuffix", &IO::PWParam::PPSuffix)
    ;

    auto c = py::class_<IO::PWPreset, IO::BasePreset>(m, "PWPreset");
    py::enum_<IO::PWPreset::AtomFmt>(c, "AtomFmt")
        .value("Bohr", IO::PWPreset::AtomFmt::Bohr)
        .value("Angstrom", IO::PWPreset::AtomFmt::Angstrom)
        .value("Crystal", IO::PWPreset::AtomFmt::Crystal)
        .value("Alat", IO::PWPreset::AtomFmt::Alat)
        .value("Active", IO::PWPreset::AtomFmt::Active)
    ;
    py::enum_<IO::PWPreset::CellFmt>(c, "CellFmt")
        .value("Angstrom", IO::PWPreset::CellFmt::Angstrom)
        .value("Bohr", IO::PWPreset::CellFmt::Bohr)
        .value("Active", IO::PWPreset::CellFmt::Active)
    ;
    c.def(py::init<IO::PWPreset::AtomFmt, IO::PWPreset::CellFmt>(),
          "atomfmt"_a=IO::PWPreset::AtomFmt::Active,
          "cellfmt"_a=IO::PWPreset::CellFmt::Active)
     .def_readwrite("atoms", &IO::PWPreset::atoms)
     .def_readwrite("cell", &IO::PWPreset::cell)
    ;
}
}
