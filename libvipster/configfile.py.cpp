#include "pyvipster.h"
#include "configfile.h"
#include "io/plugin.h"

namespace Vipster::Py{
void config(py::module& m, const ConfigState& state){
    m.attr("PeriodicTable") = std::get<0>(state);
    m.attr("Parameters") = py::cast(std::get<3>(state), py::return_value_policy::reference);
    m.attr("IOPresets") = py::cast(std::get<4>(state), py::return_value_policy::reference);
    auto plug = m.def_submodule("Plugins");
    for(const auto* p: std::get<2>(state)){
        plug.attr(p->command.c_str()) = p;
    }
}
}
