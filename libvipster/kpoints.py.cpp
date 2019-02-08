#include "pyvipster.h"
#include "kpoints.h"

namespace Vipster::Py{
void KPoints(py::module& m){
    auto k = py::class_<Vipster::KPoints>(m, "KPoints")
        .def(py::init())
        .def_readwrite("active", &KPoints::active)
        .def_readwrite("mpg", &KPoints::mpg)
        .def_readwrite("discrete", &KPoints::discrete)
    ;

    py::enum_<KPoints::Fmt>(k, "Fmt")
        .value("Gamma",KPoints::Fmt::Gamma)
        .value("MPG",KPoints::Fmt::MPG)
        .value("Discrete",KPoints::Fmt::Discrete)
    ;

    py::class_<KPoints::MPG>(k, "MPG")
        .def_readwrite("x",&KPoints::MPG::x)
        .def_readwrite("y",&KPoints::MPG::y)
        .def_readwrite("z",&KPoints::MPG::z)
        .def_readwrite("sx",&KPoints::MPG::sx)
        .def_readwrite("sy",&KPoints::MPG::sy)
        .def_readwrite("sz",&KPoints::MPG::sz)
    ;

    auto disc = py::class_<KPoints::Discrete>(k, "Discrete")
        .def_readwrite("properties", &KPoints::Discrete::properties)
        .def_readwrite("kpoints", &KPoints::Discrete::kpoints)
    ;

    py::class_<KPoints::Discrete::Point>(disc, "Point")
        .def(py::init([](const Vec& p, float w){return KPoints::Discrete::Point{p,w};}))
        .def_readwrite("pos", &KPoints::Discrete::Point::pos)
        .def_readwrite("weight", &KPoints::Discrete::Point::weight)
    ;

    py::bind_vector<std::vector<KPoints::Discrete::Point>>(disc, "Points");

    py::enum_<KPoints::Discrete::Properties>(disc, "Properties", py::arithmetic())
        .value("none", KPoints::Discrete::Properties::none)
        .value("crystal", KPoints::Discrete::Properties::crystal)
        .value("band", KPoints::Discrete::Properties::band)
        .value("contour", KPoints::Discrete::Properties::contour)
    ;

}
}
