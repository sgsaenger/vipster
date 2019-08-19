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
    ;

    py::enum_<BondPolicy>(b, "Policy")
        .value("None", BondPolicy::None)
        .value("Molecule", BondPolicy::Molecule)
        .value("Cell", BondPolicy::Cell)
    ;

    py::enum_<BondFrequency>(b, "Frequency")
        .value("Never", BondFrequency::Never)
        .value("Once", BondFrequency::Once)
        .value("Always", BondFrequency::Always)
    ;
}
}
