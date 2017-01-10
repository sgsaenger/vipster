#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <molecule.h>

int add(int i, int j) {
    return (i + 3) * j;
}

namespace py = pybind11;
using namespace py::literals;
using namespace Vipster;

PYBIND11_PLUGIN(vipster) {
    py::module m("vipster", "pybind11 example plugin");

    m.def("add", &add, "A function which adds two numbers",
            "i"_a=1, "j"_a=2);
    m.attr("what") = py::cast("Bargl?!");

    py::class_<Vec>(m, "Vec")
        .def(py::init())
        .def_property("x",[](const Vec &v){return v[0];},[](Vec &v, float f){v[0] = f;})
        .def_property("y",[](const Vec &v){return v[1];},[](Vec &v, float f){v[1] = f;})
        .def_property("z",[](const Vec &v){return v[2];},[](Vec &v, float f){v[2] = f;})
        .def("__bool__",[](const Vec&){return true;})
        //.dev("__repr__",[](....TODO
        .def(py::self == py::self);
//    py::detail::vector_if_copy_constructible; DONE
//    py::detail::vector_if_equal_operator;DONE
//    py::detail::vector_if_insertion_operator;DONE
//    py::detail::vector_modifiers; TODO setitem/getitem
//    py::detail::vector_accessor; TODO getitem/iter

    py::class_<Atom>(m, "Atom")
        .def(py::init())
        .def_readwrite("name",&Atom::name)
        .def_readwrite("coord",&Atom::coord)
        .def_readwrite("charge",&Atom::charge)
        .def_readwrite("fix",&Atom::fix)
        .def_readwrite("hidden",&Atom::hidden)
        .def(py::self == py::self);

    return m.ptr();
}
