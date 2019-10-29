#include "pyvipster.h"
#include "io/poscar/plugin.h"

namespace Vipster::Py{
void POSCAR(py::module& m){
    auto c = py::class_<IO::PoscarPreset>(m, "POSCARPreset")
        .def_readwrite("selective", &IO::PoscarPreset::selective)
        .def_readwrite("cartesian", &IO::PoscarPreset::cartesian)
    ;
}
}
