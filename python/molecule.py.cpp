#include "molecule.py.h"
#include "vipster/molecule.h"

void Vipster::Py::Molecule(py::module& m, const Vipster::ConfigState &state){
    py::class_<Vipster::Molecule>(m, "Molecule")
        .def(py::init([&state](const std::string&name, size_t s){
            Vipster::Molecule m{name, s};
            m.getPTE().root = &std::get<0>(state);
            return m;
        }), "name"_a="New Molecule", "steps"_a=1)
        .def(py::init<const Vipster::Molecule&>())
        .def(py::init<const Step&, std::string>(), "step"_a, "name"_a="Copy of Step")
        .def_property_readonly("pte", py::overload_cast<>(&Vipster::Molecule::getPTE, py::const_))
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
        .def_readwrite("name", &Vipster::Molecule::name)
        .def_readwrite("kpoints", &Vipster::Molecule::kpoints)
    ;
}
