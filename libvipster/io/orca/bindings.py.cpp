#include "pyvipster.h"
#include "io/orca/plugin.h"

namespace Vipster::Py{
void ORCA(py::module& m){
    auto p = py::class_<IO::OrcaParam, IO::BaseParam>(m, "OrcaParam")
        .def_readwrite("header", &IO::OrcaParam::header)
    ;
}
}
