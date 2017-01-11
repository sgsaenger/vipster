#include <sstream>
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <molecule.h>

namespace py = pybind11;
using namespace py::literals;
using namespace Vipster;

template <typename Array, typename holder_type = std::unique_ptr<Array>, typename... Args>
py::class_<Array, holder_type> bind_array(py::module &m, std::string const &name, Args&&... args) {
    using Class_ = pybind11::class_<Array, holder_type>;
    using ValueType = typename Array::value_type;
    using SizeType = typename Array::size_type;

    Class_ cl(m, name.c_str(), std::forward<Args>(args)...);
    cl.def(py::init());
    cl.def("__repr__",[name](const Array &v){
        std::ostringstream s;
        s << std::boolalpha << name << ": [";
        for (SizeType i=0; i < v.size()-1; ++i){
            s << v[i] << ", ";
        }
        s << v.back() << ']';
        return s.str();
    });
    cl.def("__getitem__",[](const Array &v, SizeType i){
        if(i<0 || i>= 3)
            throw py::index_error();
        return v[i];
    });
    cl.def("__setitem__",[](Array &v, SizeType i, ValueType f){
        if(i<0 || i>= 3)
            throw py::index_error();
        v[i] = f;
    });
    cl.def("__init__",[](Array &v, py::iterable it){
        new (&v) Vec();
        try{
            SizeType i=0;
            for(py::handle h : it){
                v[i] = h.cast<ValueType>();
                i++;
            }
        }catch(...){
            throw;
        }
    });
    cl.def(py::self == py::self);
    return cl;
}

PYBIND11_PLUGIN(vipster) {
    py::module m("vipster", "pybind11 example plugin");

    bind_array<Vec>(m, "Vec");
//    bind_array<Mat>(m, "Mat"); TODO: needs operator<< on c++ side...
    bind_array<FixVec>(m,"FixVec");

    py::class_<Atom>(m, "Atom")
        .def(py::init())
        .def_readwrite("name",&Atom::name)
        .def_readwrite("coord",&Atom::coord)
        .def_readwrite("charge",&Atom::charge)
        .def_readwrite("fix",&Atom::fix)
        .def_readwrite("hidden",&Atom::hidden)
        .def("__bool__",[](const Atom &a){return !a.name.empty();})
        .def(py::self == py::self);

    py::class_<Step>(m, "Step")
        .def(py::init())
        .def_property_readonly("nat", &Step::getNat)
        .def_property_readonly("ntyp", &Step::getNtyp);

    py::class_<Molecule>(m, "Molecule")
        .def(py::init());

    return m.ptr();
}
