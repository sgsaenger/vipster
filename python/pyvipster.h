#include "pybind11/pybind11.h"
#include "pybind11/operators.h"
#include "pybind11/stl_bind.h"

namespace py = pybind11;
using namespace py::literals;

template <typename Array, typename holder_type = std::unique_ptr<Array>, typename... Args>
py::class_<Array, holder_type> bind_array(py::handle &m, std::string const &name, Args&&... args) {
    using Class_ = pybind11::class_<Array, holder_type>;
    using ValueType = typename Array::value_type;
    using SizeType = typename Array::size_type;

    Class_ cl(m, name.c_str(), std::forward<Args>(args)...);
    cl.def(py::init());
    cl.def("__getitem__",[](const Array &v, SizeType i)->const ValueType&{
        if(i<0 || i>= v.size())
            throw py::index_error();
        return v[i];
    }, py::return_value_policy::reference_internal);
    cl.def("__setitem__",[](Array &v, SizeType i, ValueType f){
        if(i<0 || i>= v.size())
            throw py::index_error();
        v[i] = f;
    });
    cl.def(py::init([](const py::iterable& it){
        Array arr;
        SizeType i=0;
        for(py::handle h:it){
            arr[i] = h.cast<ValueType>();
            i++;
        }
        return arr;
    }));
    cl.def(py::self == py::self);
    py::implicitly_convertible<py::iterable, Array>();
    cl.def("__repr__", [name](const Array &v){
        std::string repr = name + "[";
        for(size_t i=0; i<v.size()-1; ++i){
            repr += py::str(py::cast(v[i]).attr("__repr__")());
            repr += ", ";
        }
        repr += py::str(py::cast(v[v.size()-1]).attr("__repr__")());
        repr += "]";
        return repr;
    });
    return cl;
}
