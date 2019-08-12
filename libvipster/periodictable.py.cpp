#include "pyvipster.h"
#include "periodictable.h"

namespace Vipster::Py{
void Table(py::module& m){
    auto p = py::class_<PeriodicTable, std::shared_ptr<PeriodicTable>>(m, "PeriodicTable")
        .def("__getitem__", &PeriodicTable::operator [], py::return_value_policy::reference_internal)
        .def("__setitem__", [](PeriodicTable &pte, std::string n, Element& e){pte[n] = e;})
    ;

    py::class_<Element>(p, "Element")
        .def_readwrite("PWPP", &Element::PWPP)
        .def_readwrite("CPPP", &Element::CPPP)
        .def_readwrite("CPNL", &Element::CPNL)
        .def_readwrite("Z", &Element::Z)
        .def_readwrite("m", &Element::m)
        .def_readwrite("bondcut", &Element::bondcut)
        .def_readwrite("covr", &Element::covr)
        .def_readwrite("vdwr", &Element::vdwr)
        .def_readwrite("col", &Element::col)
    ;
}
}
