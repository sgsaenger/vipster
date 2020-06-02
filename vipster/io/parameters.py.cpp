#include "parameters.py.h"
#include "plugin.h"
#include "parameters.h"

namespace Vipster::IO {
std::ostream& operator<<(std::ostream& os, const Parameter &p){
    nlohmann::json j;
    to_json(j, p);
    return os << j;
}
}

PYBIND11_MAKE_OPAQUE(Vipster::IO::Parameters);
PYBIND11_MAKE_OPAQUE(Vipster::IO::Parameters::mapped_type);

void Vipster::Py::Parameters(py::module &io){
    py::class_<IO::Parameter>(io, "Parameter")
        .def("__repr__", [](const IO::Parameter &p){
            std::ostringstream s;
            nlohmann::json j;
            to_json(j, p);
            s << p.getFmt()->command << "-Parameter" << j;
            return s.str();
        });

    py::bind_map<IO::Parameters::mapped_type>(io, "__StrParMap__");
    py::bind_map<IO::Parameters>(io, "Parameters")
        .def("__repr__", [](IO::Parameters& p){
            std::ostringstream s;
            s << "Parameters{";
            bool f = false;
            for (auto const &kv:p){
                if (f)
                    s << ", ";
                s << kv.first << ": __StrParMap__";
                f = true;
            }
            s << '}';
            return s.str();
        })
    ;
}
