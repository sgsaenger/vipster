#include "pybind11/pybind11.h"
#include "pybind11/operators.h"
#include "pybind11/stl_bind.h"

#include <sstream>
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
    return cl;
}

PYBIND11_MODULE(vipster, m) {
    m.doc() = "Vipster python bindings";

    /*
     * Step (Atoms, Bonds, Cell + PSE)
     */

    bind_array<Vec>(m, "Vec");
    bind_array<Mat>(m, "Mat");
    bind_array<ColVec>(m, "ColVec");
    py::bind_vector<std::vector<StepProper>>(m,"__StepVector__");
    py::bind_vector<std::vector<Bond>>(m,"__BondVector__");
    py::bind_vector<std::vector<std::string>>(m,"__StringVector__");

    py::enum_<AtomFmt>(m, "AtomFmt")
        .value("Bohr", AtomFmt::Bohr)
        .value("Angstrom", AtomFmt::Angstrom)
        .value("Crystal", AtomFmt::Crystal)
        .value("Alat", AtomFmt::Alat)
    ;

    py::enum_<BondLevel>(m, "BondLevel")
        .value("None", BondLevel::None)
        .value("Molecule", BondLevel::Molecule)
        .value("Cell", BondLevel::Cell)
    ;

    py::enum_<CdmFmt>(m, "CellDimFmt")
        .value("Bohr", CdmFmt::Bohr)
        .value("Angstrom", CdmFmt::Angstrom)
    ;

    py::class_<Atom>(m, "Atom")
        .def_property("name", [](const Atom &a)->const std::string&{return a.name;},
                      [](Atom &a, std::string s){a.name = s;})
        .def_property("coord", [](const Atom &a)->const Vec&{return a.coord;},
                      [](Atom &a, Vec c){a.coord = c;})
        .def_property("charge", [](const Atom &a)->const float&{return a.charge;},
                      [](Atom &a, float c){a.charge = c;})
    //TODO: make bitset observable. Should pse-entry be exposed?
//        .def_property("properties", [](const Atom &a)->const )
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
        .def_readonly("pse", &Step::pse)
        .def_property("comment", &Step::getComment, &Step::setComment)
        .def("newAtom", [](Step& s){s.newAtom();})
    //TODO: enable when bitset is wrapped
//        .def("newAtom", py::overload_cast<std::string, Vec, float, std::bitset<nAtProp>>(&StepProper::newAtom),
//             "name"_a, "coord"_a=Vec{}, "charge"_a=float{})
        .def("newAtom", py::overload_cast<const Atom&>(&StepProper::newAtom), "at"_a)
        .def("__getitem__", [](Step& s, int i){
                if (i<0){
                    i = i+static_cast<int>(s.getNat());
                }
                if ((i<0) || i>=static_cast<int>(s.getNat())){
                    throw py::index_error();
                }
                return s[static_cast<size_t>(i)];
            })
        .def("__setitem__", [](Step& s, int i, const Atom& at){
                if (i<0){
                    i = i+static_cast<int>(s.getNat());
                }
                if ((i<0) || i>=static_cast<int>(s.getNat())){
                    throw py::index_error();
                }
                s[static_cast<size_t>(i)] = at;
        })
        .def("__iter__", [](const Step& s){return py::make_iterator(s.begin(), s.end());})
        .def("delAtom", &Step::delAtom, "i"_a)
        .def_property_readonly("nat", &Step::getNat)
        .def("getTypes", [](const Step& s){
            auto oldT = s.getTypes();
            return std::vector<std::string>(oldT.begin(), oldT.end());
        })
        .def_property_readonly("ntyp", &Step::getNtyp)
    //FMT
        .def("getFmt", &Step::getFmt)
        .def("asFmt", py::overload_cast<AtomFmt>(&Step::asFmt), "fmt"_a,
             py::return_value_policy::reference_internal)
    //CELL
        .def("enableCell", &Step::enableCell, "enable"_a)
        .def("getCellDim", &Step::getCellDim)
        .def("setCellDim", &Step::setCellDim, "cdm"_a, "fmt"_a, "scale"_a=false)
        .def("getCellVec", &Step::getCellVec)
        .def("setCellVec", &Step::setCellVec, "vec"_a, "scale"_a=false)
        .def("getCenter", &Step::getCenter, "fmt"_a, "com"_a=false)
    //BONDS
        .def("getBonds", py::overload_cast<BondLevel, bool>(&Step::getBonds, py::const_),
             "level"_a=BondLevel::Cell, "update"_a=true)
        .def("getBonds", py::overload_cast<float, BondLevel, bool>(&StepProper::getBonds, py::const_),
             "cutfac"_a, "level"_a=BondLevel::Cell, "update"_a=true)
        .def_property_readonly("nbond", &StepProper::getNbond)
    ;

    //TODO: allow construction?
    py::class_<StepProper, Step>(m, "StepProper")
        .def("setFmt", &StepProper::setFmt, "fmt"_a, "scale"_a=false);
    py::class_<StepFormatter, Step>(m, "StepFormatter");

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

    /*
     * Molecule (Steps + KPoints)
     */
    py::class_<Molecule>(m, "Molecule")
        .def(py::init())
        .def_readonly("pse", &Molecule::pse)
        .def("newStep", [](Molecule& m){m.newStep();})
        .def("newStep", py::overload_cast<const StepProper&>(&Molecule::newStep), "step"_a)
        .def("getStep", py::overload_cast<size_t>(&Molecule::getStep), "i"_a, py::return_value_policy::reference_internal)
        .def("getSteps",py::overload_cast<>(&Molecule::getSteps))
        .def_property_readonly("nstep", &Molecule::getNstep)
        .def_property("name", &Molecule::getName, &Molecule::setName)
        .def_property("kpoints", py::overload_cast<>(&Molecule::getKPoints), &Molecule::setKPoints)
    ;

    /*
     * IO
     */
    py::enum_<IOFmt>(m,"IOFmt")
        .value("XYZ", IOFmt::XYZ)
        .value("PWI", IOFmt::PWI)
        .value("PWO", IOFmt::PWO)
        .value("LMP", IOFmt::LMP)
        .value("DMP", IOFmt::DMP)
    ;

    m.def("readFile",[](std::string fn, IOFmt fmt){return readFile(fn,fmt).mol;},"fn"_a,"fmt"_a);
}
