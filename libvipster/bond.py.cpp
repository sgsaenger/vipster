#include "pyvipster.h"
#include "bond.h"

namespace Vipster::Py{
void Bond(py::module &m){
    py::bind_vector<std::vector<Vipster::Bond>>(m,"__BondVector__");

    auto b = py::class_<Vipster::Bond>(m, "Bond")
        .def_readwrite("at1", &Bond::at1)
        .def_readwrite("at2", &Bond::at2)
        .def_readwrite("dist", &Bond::dist)
        .def_readwrite("diff", &Bond::diff)
        .def_property_readonly("type", [](const Vipster::Bond &b){return b.type->first;})
    ;

    py::enum_<BondMode>(m, "Mode")
        .value("Manual", BondMode::Manual)
        .value("Automatic", BondMode::Automatic)
    ;
}
}
