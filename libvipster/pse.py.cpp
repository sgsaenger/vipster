#include "pyvipster.h"
#include "pse.h"

namespace Vipster::Py{
void PSE(py::module& m){
    auto p = py::class_<PseMap, std::shared_ptr<PseMap>>(m, "PseMap")
        .def("__getitem__", &PseMap::operator [], py::return_value_policy::reference_internal)
        .def("__setitem__", [](PseMap &pse, std::string n, PseEntry& e){pse[n] = e;})
    ;

    py::class_<PseEntry>(p, "Entry")
        .def_readwrite("PWPP", &PseEntry::PWPP)
        .def_readwrite("CPPP", &PseEntry::CPPP)
        .def_readwrite("CPNL", &PseEntry::CPNL)
        .def_readwrite("Z", &PseEntry::Z)
        .def_readwrite("m", &PseEntry::m)
        .def_readwrite("bondcut", &PseEntry::bondcut)
        .def_readwrite("covr", &PseEntry::covr)
        .def_readwrite("vdwr", &PseEntry::vdwr)
        .def_readwrite("col", &PseEntry::col)
    ;
}
}
