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
        s << v;
        return s.str();
    });
    cl.def("__getitem__",[](const Array &v, SizeType i){
        if(i<0 || i>= v.size())
            throw py::index_error();
        return v[i];
    });
    cl.def("__setitem__",[](Array &v, SizeType i, ValueType f){
        if(i<0 || i>= v.size())
            throw py::index_error();
        v[i] = f;
    });
    cl.def("__init__",[](Array &v, py::iterable it){
        new (&v) Array();
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
    bind_array<Mat>(m, "Mat");
    bind_array<FixVec>(m,"FixVec");
    py::implicitly_convertible<py::iterable, Vec>();

    py::enum_<AtomFmt>(m, "AtomFmt")
        .value("Bohr", AtomFmt::Bohr)
        .value("Angstrom", AtomFmt::Angstrom)
        .value("Crystal", AtomFmt::Crystal)
        .value("Alat", AtomFmt::Alat);

    py::class_<Atom>(m, "Atom")
        .def("__init__",
             [](Atom &a, std::string n, Vec coord, float charge, FixVec fix, bool hidden)
             { new (&a) Atom{n,coord,charge,fix,hidden};  },
             "name"_a="C", "coord"_a=Vec{{0,0,0}}, "charge"_a=0,
             "fix"_a=FixVec{{false,false,false}}, "hidden"_a=false )
        .def_readwrite("name",&Atom::name)
        .def_readwrite("coord",&Atom::coord)
        .def_readwrite("charge",&Atom::charge)
        .def_readwrite("fix",&Atom::fix)
        .def_readwrite("hidden",&Atom::hidden)
        .def("__bool__",[](const Atom &a){return !a.name.empty();})
        .def("__repr__",[](const Atom &a){std::ostringstream s;s<<a;return s.str();})
        .def(py::self == py::self)
        .def(py::self != py::self)
    ;

    py::class_<Step>(m, "Step")
        .def(py::init())
        .def("__repr__",[](const Step &a){std::ostringstream s;s<<a;return s.str();})
        .def_property_readonly("nat", &Step::getNat)
        .def_property_readonly("ntyp", &Step::getNtyp)
//        .def("newAtom", [](Step& s){s.newAtom();})
//        .def("newAtom", py::overload_cast<const Atom&>(&Step::newAtom), "at"_a)
        .def("newAtom", py::overload_cast<Atom, AtomFmt>(&Step::newAtom), "at"_a=Atom{"C"}, "fmt"_a=AtomFmt::Bohr)
        .def("delAtom", &Step::delAtom, "i"_a)
        .def("setAtom", py::overload_cast<size_t, const Atom&>(&Step::setAtom), "i"_a, "at"_a)
        .def("setAtom", py::overload_cast<size_t, Atom, AtomFmt>(&Step::setAtom), "i"_a, "at"_a, "fmt"_a)
        .def("getAtom", &Step::getAtom, "i"_a, "fmt"_a=AtomFmt::Bohr)
        .def("getCellDim", &Step::getCellDim, "fmt"_a=AtomFmt::Bohr)
        .def("setCellDim", &Step::setCellDim, "dim"_a, "scale"_a=false, "fmt"_a=AtomFmt::Bohr)
        .def("getCenter", &Step::getCenter, "com"_a=false)
        .def_property("comment", &Step::getComment, &Step::setComment)
    ;

    py::class_<Molecule>(m, "Molecule")
        .def(py::init())
        .def("__repr__",[](const Molecule &m){std::ostringstream s;s<<m;return s.str();})
    ;

    return m.ptr();
}
