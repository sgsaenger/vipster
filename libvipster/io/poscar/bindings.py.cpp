#include "pyvipster.h"
#include "io/poscar/plugin.h"

namespace Vipster::Py{
void POSCAR(py::module& m){
    auto c = py::class_<IO::PoscarPreset, IO::BasePreset>(m, "POSCARPreset")
        .def(py::init<bool, bool>(), "selective"_a=true, "cartesian"_a=false)
        .def_readwrite("selective", &IO::PoscarPreset::selective)
        .def_readwrite("cartesian", &IO::PoscarPreset::cartesian)
    ;
}
}
