#include <sstream>
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
//#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <molecule.h>
#include <iowrapper.h>

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
    py::implicitly_convertible<py::iterable, Array>();
    return cl;
}

PYBIND11_PLUGIN(vipster) {
    py::module m("vipster", "pybind11 example plugin");

    /*
     * Step (Atoms, Bonds, Cell + PSE)
     */

    bind_array<Vec>(m, "Vec");
    bind_array<Mat>(m, "Mat");
    bind_array<FixVec>(m, "FixVec");
    bind_array<ColVec>(m, "ColVec");
    py::bind_vector<std::vector<Atom>>(m,"__AtomVector__");
    py::bind_vector<std::vector<Step>>(m,"__StepVector__");
    py::bind_vector<std::vector<Bond>>(m,"__BondVector__");

    py::enum_<AtomFmt>(m, "AtomFmt")
        .value("Bohr", AtomFmt::Bohr)
        .value("Angstrom", AtomFmt::Angstrom)
        .value("Crystal", AtomFmt::Crystal)
        .value("Alat", AtomFmt::Alat)
    ;

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
        .def(py::self == py::self)
        .def(py::self != py::self)
    ;

    py::class_<Bond>(m, "Bond")
        .def_readwrite("at1", &Bond::at1)
        .def_readwrite("at2", &Bond::at2)
        .def_readwrite("dist", &Bond::dist)
        .def_readwrite("xdiff", &Bond::xdiff)
        .def_readwrite("ydiff", &Bond::ydiff)
        .def_readwrite("zdiff", &Bond::zdiff)
    ;

    py::class_<PseMap, std::shared_ptr<PseMap>>(m, "PseMap")
        .def("__getitem__", &PseMap::operator [])
        .def("__setitem__", [](PseMap &pse, std::string n, PseEntry& e){pse[n] = e;})
    ;

    py::class_<PseEntry>(m, "PseEntry")
        .def_readwrite("PWPP", &PseEntry::PWPP)
        .def_readwrite("CPPP", &PseEntry::CPPP)
        .def_readwrite("CPNL", &PseEntry::CPNL)
        .def_readwrite("Z", &PseEntry::Z)
        .def_readwrite("m", &PseEntry::m)
        .def_readwrite("bondcut", &PseEntry::bondcut)
        .def_readwrite("covr", &PseEntry::covr)
        .def_readwrite("vdwr", &PseEntry::vdwr)
        .def_readwrite("col", &PseEntry::col)
    ;

    py::class_<Step>(m, "Step")
        .def(py::init())
        .def_readonly("pse", &Step::pse)
        .def_readwrite("comment", &Step::comment)
        .def("newAtom", [](Step& s){s.newAtom();})
        .def("newAtom", py::overload_cast<const Atom&>(&Step::newAtom), "at"_a)
        .def("newAtom", py::overload_cast<Atom, AtomFmt>(&Step::newAtom), "at"_a, "fmt"_a)
        .def("delAtom", &Step::delAtom, "i"_a)
        .def("setAtom", py::overload_cast<size_t, const Atom&>(&Step::setAtom), "i"_a, "at"_a)
        .def("setAtom", py::overload_cast<size_t, Atom, AtomFmt>(&Step::setAtom), "i"_a, "at"_a, "fmt"_a)
        .def("getAtom", py::overload_cast<size_t>(&Step::getAtom, py::const_), "i"_a)
        .def("getAtomFmt", py::overload_cast<size_t, AtomFmt>(&Step::getAtomFmt, py::const_), "i"_a, "fmt"_a)
        .def("getAtoms", py::overload_cast<>(&Step::getAtoms, py::const_))
        .def("getAtomsFmt", py::overload_cast<AtomFmt>(&Step::getAtomsFmt, py::const_), "fmt"_a)
        .def_property_readonly("nat", &Step::getNat)
        .def("getTypes", &Step::getTypes)
        .def_property_readonly("ntyp", &Step::getNtyp)
        .def("getFmt", &Step::getFmt)
        .def("setFmt", &Step::setFmt, "fmt"_a, "scale"_a=false)
        .def("getCellDim", py::overload_cast<>(&Step::getCellDim, py::const_))
        .def("getCellDim", py::overload_cast<AtomFmt>(&Step::getCellDim, py::const_), "fmt"_a)
        .def("setCellDim", py::overload_cast<float,bool>(&Step::setCellDim), "dim"_a, "scale"_a=false)
        .def("setCellDim", py::overload_cast<float,bool,AtomFmt>(&Step::setCellDim), "dim"_a, "scale"_a, "fmt"_a)
        .def("getCellVec", &Step::getCellVec)
        .def("setCellVec", &Step::setCellVec, "vec"_a, "scale"_a=false)
        .def("getCenter", &Step::getCenter, "com"_a=false)
        .def("getBonds", py::overload_cast<>(&Step::getBonds, py::const_))
        .def("getBonds", py::overload_cast<float>(&Step::getBonds, py::const_), "cutfac"_a)
        .def("getBondsCell", py::overload_cast<>(&Step::getBondsCell, py::const_))
        .def("getBondsCell", py::overload_cast<float>(&Step::getBondsCell, py::const_), "cutfac"_a)
        .def_property_readonly("nbond", &Step::getNbond)
    ;

    /*
     * K-Points
     */
    py::bind_vector<std::vector<DiscreteKPoint>>(m, "__KPointVector__");

    py::enum_<KPointFmt>(m, "KPointFmt")
        .value("Gamma",KPointFmt::Gamma)
        .value("MPG",KPointFmt::MPG)
        .value("Discrete",KPointFmt::Discrete)
    ;

    py::class_<KPoints> k(m, "KPoints");
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
    dp.def("__init__",[](DiscreteKPoint& dp, const Vec& p, float w)
           { new (&dp) DiscreteKPoint{p,w};},
           "pos"_a=Vec{},"weight"_a=0.)
        .def_readwrite("pos", &DiscreteKPoint::pos)
        .def_readwrite("weight", &DiscreteKPoint::weight)
    ;

    py::enum_<KPoints::Discrete::Properties>(dp, "Properties", py::arithmetic())
        .value("none", KPoints::Discrete::Properties::none)
        .value("crystal", KPoints::Discrete::Properties::crystal)
        .value("band", KPoints::Discrete::Properties::band)
    ;

    /*
     * Molecule (Steps + KPoints)
     */
    py::class_<Molecule>(m, "Molecule")
        .def(py::init())
        .def("newStep", py::overload_cast<const Step&>(&Molecule::newStep), "step"_a=Step{})
        .def("getStep", py::overload_cast<size_t>(&Molecule::getStep), "i"_a)
        .def("getSteps",py::overload_cast<>(&Molecule::getSteps))
        .def_property_readonly("nstep", &Molecule::getNstep)
        .def_property("name", &Molecule::getName, &Molecule::setName)
        .def_property("kpoints", py::overload_cast<>(&Molecule::getKPoints), &Molecule::setKPoints)
    ;

    /*
     * IO
     */
    py::enum_<IOFmt>(m,"IOFmt")
        .value("XYZ",IOFmt::XYZ)
        .value("PWI",IOFmt::PWI)
    ;

    m.def("readFile",[](std::string fn, IOFmt fmt){return readFile(fn,fmt)->mol;},"fn"_a,"fmt"_a);

    return m.ptr();
}
