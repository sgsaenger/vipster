#include "pyvipster.h"
#include "vec.h"

namespace Vipster::Py{
void Vec(py::module& m){
    bind_array<Vipster::Vec>(m, "Vec")
        // TODO: shorthand definitions?
        .def("__add__", [](const Vipster::Vec &v1, const Vipster::Vec &v2){return v1+v2;}, py::is_operator())
        .def("__sub__", [](const Vipster::Vec &v1, const Vipster::Vec &v2){return v1-v2;}, py::is_operator())
        .def("__add__", [](const Vipster::Vec &v, float f){return v+f;}, py::is_operator())
        .def("__sub__", [](const Vipster::Vec &v, float f){return v-f;}, py::is_operator())
        .def("__mul__", [](const Vipster::Vec &v, float f){return v*f;}, py::is_operator())
        .def("__truediv__", [](const Vipster::Vec &v, float f){return v/f;}, py::is_operator())
        .def("__radd__", [](const Vipster::Vec &v, float f){return f+v;}, py::is_operator())
        .def("__rmul__", [](const Vipster::Vec &v, float f){return f*v;}, py::is_operator())
        .def("dot", &Vec_dot)
        .def("cross", &Vec_cross)
        .def("len", &Vec_length)
    ;
    bind_array<Mat>(m, "Mat")
        .def("__mul__", [](const Mat &m, float f){return m*f;}, py::is_operator())
        .def("__rmul__", [](const Mat &m, float f){return f*m;}, py::is_operator())
        .def("__truediv__", [](const Mat &m, float f){return m/f;}, py::is_operator())
        .def("__mul__", [](const Mat &m, const Vipster::Vec& v){return m*v;}, py::is_operator())
        .def("__rmul__", [](const Mat &m, const Vipster::Vec& v){return v*m;}, py::is_operator())
        .def("__mul__", [](const Mat &m1, const Mat &m2){return m1*m2;}, py::is_operator())
        .def("transpose", &Mat_trans)
        .def("determinant", &Mat_det)
        .def("inverse", &Mat_inv)
    ;
}
}
