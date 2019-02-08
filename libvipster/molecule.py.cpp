#include "pyvipster.h"
#include "molecule.h"

namespace Vipster::Py{
void Molecule(py::module& m){
    py::class_<Vipster::Molecule>(m, "Molecule")
        .def(py::init())
        .def_readonly("pse", &Vipster::Molecule::pse)
        .def("newStep", [](Vipster::Molecule& m){m.newStep();})
        .def("newStep", py::overload_cast<const Step&>(&Vipster::Molecule::newStep), "step"_a)
        .def("__getitem__", [](Vipster::Molecule& m, int i)->Step&{
            if(i<0){
                i = i+static_cast<int>(m.getNstep());
            }
            if((i<0) || i>=static_cast<int>(m.getNstep())){
                throw py::index_error();
            }
            return m.getStep(static_cast<size_t>(i));
        }, py::return_value_policy::reference_internal)
        .def("__iter__", [](const Vipster::Molecule& m){return py::make_iterator(m.getSteps().begin(), m.getSteps().end());})
        .def("__len__", &Vipster::Molecule::getNstep)
        .def_property_readonly("nstep", &Vipster::Molecule::getNstep)
        .def_property("name", &Vipster::Molecule::getName, &Vipster::Molecule::setName)
        .def_property("kpoints", py::overload_cast<>(&Vipster::Molecule::getKPoints), &Vipster::Molecule::setKPoints)
    ;
}
}
