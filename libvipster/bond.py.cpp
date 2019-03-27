#include "pyvipster.h"
#include "bond.h"

namespace Vipster::Py{
void Bond(py::module &m){
    py::bind_vector<std::vector<Vipster::Bond>>(m,"__BondVector__");

    auto b = py::class_<Vipster::Bond>(m, "Bond")
        .def_readwrite("at1", &Bond::at1)
        .def_readwrite("at2", &Bond::at2)
        .def_readwrite("dist", &Bond::dist)
        .def_readwrite("xdiff", &Bond::xdiff)
        .def_readwrite("ydiff", &Bond::ydiff)
        .def_readwrite("zdiff", &Bond::zdiff)
    ;

    py::enum_<BondLevel>(b, "Level")
        .value("None", BondLevel::None)
        .value("Molecule", BondLevel::Molecule)
        .value("Cell", BondLevel::Cell)
    ;

    py::enum_<BondFrequency>(b, "Frequency")
        .value("Never", BondFrequency::Never)
        .value("Once", BondFrequency::Once)
        .value("Always", BondFrequency::Always)
    ;
}
}
