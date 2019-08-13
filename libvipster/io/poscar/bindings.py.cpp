#include "pyvipster.h"
#include "io/poscar/plugin.h"

namespace Vipster::Py{
void POSCAR(py::module& m){
    auto c = py::class_<IO::PoscarConfig>(m, "POSCARConfig")
        .def_readwrite("selective", &IO::PoscarConfig::selective)
        .def_readwrite("cartesian", &IO::PoscarConfig::cartesian)
    ;
}
}
