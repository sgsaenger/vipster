#include "pyvipster.h"
#include "bond.h"
#include <pybind11/stl.h>

namespace Vipster{
std::ostream& operator<<(std::ostream& os, const Vipster::Bond&b){
    os << "Bond(" << b.at1 << ", " << b.at2 << ", <"
       << b.diff[0] << ", " << b.diff[1] << ", " << b.diff[2]
       << '>';
    if(b.type){
        os << ", " << b.type->first;
    }
    os << ')';
    return os;
}
}

PYBIND11_MAKE_OPAQUE(std::vector<Vipster::Bond>);

namespace Vipster::Py{
void Bond(py::module &m){
    py::bind_vector<std::vector<Vipster::Bond>>(m,"__BondVector__");

    auto b = py::class_<Vipster::Bond>(m, "Bond")
        .def("__repr__", [](const Vipster::Bond&b){
            std::ostringstream s;
            s << b;
            return s.str();
        })
        .def_readwrite("at1", &Bond::at1)
        .def_readwrite("at2", &Bond::at2)
        .def_readwrite("dist", &Bond::dist)
        .def_readwrite("diff", &Bond::diff)
        .def_property_readonly("type", [](const Vipster::Bond &b)->std::string{
            if(b.type){
                return b.type->first;
            }else{
                return "";
            }
        })
    ;

    py::enum_<BondMode>(m, "Mode")
        .value("Manual", BondMode::Manual)
        .value("Automatic", BondMode::Automatic)
    ;
}
}
