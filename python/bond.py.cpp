#include "bond.py.h"
#include "vipster/bond.h"
#include <pybind11/stl.h>

PYBIND11_MAKE_OPAQUE(std::vector<Vipster::Bond>)
PYBIND11_MAKE_OPAQUE(std::vector<Vipster::Overlap>)
PYBIND11_MAKE_OPAQUE(std::vector<Vipster::Angle>)
PYBIND11_MAKE_OPAQUE(std::vector<Vipster::Dihedral>)

void Vipster::Py::Bond(py::module &m){
    py::bind_vector<std::vector<Vipster::Bond>>(m,"__BondVector__");
    py::bind_vector<std::vector<Overlap>>(m,"__OverlapVector__");
    py::bind_vector<std::vector<Angle>>(m,"__AngleVector__");
    py::bind_vector<std::vector<Dihedral>>(m,"__DihedralVector__");

    py::class_<Vipster::Bond>(m, "Bond")
        .def("__repr__", [](const Vipster::Bond&b){
            std::ostringstream os;
            os << "Bond(" << b.at1 << ", " << b.at2 << ", <"
               << b.diff[0] << ", " << b.diff[1] << ", " << b.diff[2]
               << '>';
            if(b.type){
                os << ", " << b.type->first;
            }
            os << ')';
            return os.str();
        })
        .def_readonly("at1", &Bond::at1)
        .def_readonly("at2", &Bond::at2)
        .def_readonly("dist", &Bond::dist)
        .def_readonly("diff", &Bond::diff)
        .def_property_readonly("type", [](const Vipster::Bond &b)->std::string{
            if(b.type){
                return b.type->first;
            }else{
                return "";
            }
        })
    ;

    py::class_<Overlap>(m, "Overlap")
        .def("__repr__", [](const Overlap &o){
            std::ostringstream os;
            os << "Overlap(" << o.at1 << ", " << o.at2 << ')';
            return os.str();
        })
        .def_readonly("at1", &Overlap::at1)
        .def_readonly("at2", &Overlap::at2)
        .def_readonly("periodic", &Overlap::periodic)
    ;

    py::class_<Angle>(m, "Angle")
        .def("__repr__", [](const Angle &a){
            std::ostringstream os;
            os << "Angle(" << a.at1 << ", " << a.at2 << ", " << a.at3 << ')';
            return os.str();
        })
        .def_readonly("at1", &Angle::at1)
        .def_readonly("at2", &Angle::at2)
        .def_readonly("at3", &Angle::at3)
    ;

    py::class_<Dihedral>(m, "Dihedral")
        .def("__repr__", [](const Dihedral &a){
            std::ostringstream os;
            os << "Dihedral(" << a.at1 << ", " << a.at2 << ", " << a.at3 << ')';
            return os.str();
        })
        .def_readonly("at1", &Dihedral::at1)
        .def_readonly("at2", &Dihedral::at2)
        .def_readonly("at3", &Dihedral::at3)
        .def_readonly("at4", &Dihedral::at4)
    ;
}
