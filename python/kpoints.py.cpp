#include "pyvipster.h"
#include "kpoints.h"

namespace Vipster{
namespace Py{
void KPoints(py::module& m){
    py::bind_vector<std::vector<DiscreteKPoint>>(m, "__KPointVector__");

    py::enum_<KPointFmt>(m, "KPointFmt")
        .value("Gamma",KPointFmt::Gamma)
        .value("MPG",KPointFmt::MPG)
        .value("Discrete",KPointFmt::Discrete)
    ;

    py::class_<Vipster::KPoints> k(m, "KPoints");
    k.def(py::init())
        .def_readwrite("active", &KPoints::active)
        .def_readwrite("mpg", &KPoints::mpg)
        .def_readwrite("discrete", &KPoints::discrete)
    ;

    py::class_<KPoints::MPG>(k, "MPG")
        .def_readwrite("x",&KPoints::MPG::x)
        .def_readwrite("y",&KPoints::MPG::y)
        .def_readwrite("z",&KPoints::MPG::z)
        .def_readwrite("sx",&KPoints::MPG::sx)
        .def_readwrite("sy",&KPoints::MPG::sy)
        .def_readwrite("sz",&KPoints::MPG::sz)
    ;

    py::class_<KPoints::Discrete> disc(k, "Discrete");
    disc.def_readwrite("properties", &KPoints::Discrete::properties)
        .def_readwrite("kpoints", &KPoints::Discrete::kpoints)
    ;

    py::class_<DiscreteKPoint> dp(m, "DiscreteKPoint");
    dp.def(py::init([](const Vec& p, float w){return DiscreteKPoint{p,w};}))
        .def_readwrite("pos", &DiscreteKPoint::pos)
        .def_readwrite("weight", &DiscreteKPoint::weight)
    ;

    py::enum_<KPoints::Discrete::Properties>(dp, "Properties", py::arithmetic())
        .value("none", KPoints::Discrete::Properties::none)
        .value("crystal", KPoints::Discrete::Properties::crystal)
        .value("band", KPoints::Discrete::Properties::band)
        .value("contour", KPoints::Discrete::Properties::contour)
    ;

}
}
}
