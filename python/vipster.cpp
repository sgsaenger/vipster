#include <pybind11/pybind11.h>
#include <vipster/molecule.h>

int add(int i, int j) {
    return (i + 3) * j;
}

namespace py = pybind11;
using namespace py::literals;

PYBIND11_PLUGIN(test) {
    py::module m("test", "pybind11 example plugin");

    m.def("add", &add, "A function which adds two numbers",
            "i"_a=1, "j"_a=2);
    m.attr("what") = py::cast("Bargl?!");

    py::class_<Vipster::Atom>(m, "Atom");

    return m.ptr();
}
