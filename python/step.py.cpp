#include "pyvipster.h"
#include "step.h"

// TODO: think about const-correctness!
// can wait until "view"-interface to running libvipster instance is implemented

namespace Vipster{
namespace Py{
void Step(py::module& m){
    py::enum_<CdmFmt>(m, "CellDimFmt")
        .value("Bohr", CdmFmt::Bohr)
        .value("Angstrom", CdmFmt::Angstrom)
    ;

    py::class_<Vipster::Step>(m, "Vipster::Step")
        .def(py::init<AtomFmt, std::string>(), "fmt"_a=AtomFmt::Bohr, "comment"_a="")
        .def_readonly("pse", &Vipster::Step::pse)
        .def_property("comment", &Vipster::Step::getComment, &Vipster::Step::setComment)
    // Atoms
        .def("newAtom", py::overload_cast<std::string, Vec, AtomProperties>(&Vipster::Step::newAtom),
             "name"_a="", "coord"_a=Vec{}, "properties"_a=AtomProperties{})
        .def("newAtom", py::overload_cast<const Atom&>(&Vipster::Step::newAtom), "at"_a)
        .def("__getitem__", [](Vipster::Step& s, int i){
                if (i<0){
                    i = i+static_cast<int>(s.getNat());
                }
                if ((i<0) || i>=static_cast<int>(s.getNat())){
                    throw py::index_error();
                }
                return s[static_cast<size_t>(i)];
            })
        .def("__setitem__", [](Vipster::Step& s, int i, const Atom& at){
                if (i<0){
                    i = i+static_cast<int>(s.getNat());
                }
                if ((i<0) || i>=static_cast<int>(s.getNat())){
                    throw py::index_error();
                }
                s[static_cast<size_t>(i)] = at;
        })
        .def("__len__", &Vipster::Step::getNat)
        .def_property_readonly("nat", &Vipster::Step::getNat)
        .def("__iter__", [](const Vipster::Step& s){return py::make_iterator(s.begin(), s.end());})
        .def("delAtom", &Vipster::Step::delAtom, "i"_a)
    // TYPES
        .def("getTypes", [](const Vipster::Step& s){
            auto oldT = s.getTypes();
            return std::vector<std::string>(oldT.begin(), oldT.end());
        })
        .def_property_readonly("ntyp", &Vipster::Step::getNtyp)
    // FMT
        .def("getFmt", &Vipster::Step::getFmt)
        .def("asFmt", py::overload_cast<AtomFmt>(&Vipster::Step::asFmt), "fmt"_a,
             py::return_value_policy::reference_internal)
        .def("setFmt", &Vipster::Step::setFmt, "fmt"_a)
    // CELL
        .def("hasCell", &Vipster::Step::hasCell)
        .def("enableCell", &Vipster::Step::enableCell, "enable"_a)
        .def("getCellDim", &Vipster::Step::getCellDim)
        .def("setCellDim", &Vipster::Step::setCellDim, "cdm"_a, "fmt"_a, "scale"_a=false)
        .def("getCellVec", &Vipster::Step::getCellVec)
        .def("setCellVec", &Vipster::Step::setCellVec, "vec"_a, "scale"_a=false)
        .def("getCom", py::overload_cast<>(&Vipster::Step::getCom, py::const_))
        .def("getCom", py::overload_cast<AtomFmt>(&Vipster::Step::getCom, py::const_), "fmt"_a)
        .def("getCenter", &Vipster::Step::getCenter, "fmt"_a, "com"_a=false)
    // BONDS
        .def("getBonds", &Vipster::Step::getBonds,
             "cutfac"_a=settings.bondCutFac.val, "level"_a=settings.bondLvl.val, "update"_a=settings.bondFreq.val)
        .def("setBonds", &Vipster::Step::setBonds,
             "level"_a=settings.bondLvl.val, "cutfac"_a=settings.bondCutFac.val)
        .def_property_readonly("nbond", &Vipster::Step::getNbond)
    // SELECTION
        .def("select", py::overload_cast<std::string>(&Vipster::Step::select), "filter"_a)
    // TODO: Modification functions
    ;

    py::class_<Step::selection>(m, "Selection")
        .def_readonly("pse", &Step::selection::pse)
        .def_property("comment", &Step::selection::getComment, &Step::selection::setComment)
    // Atoms
        .def("__getitem__", [](Step::selection& s, int i){
                if (i<0){
                    i = i+static_cast<int>(s.getNat());
                }
                if ((i<0) || i>=static_cast<int>(s.getNat())){
                    throw py::index_error();
                }
                return s[static_cast<size_t>(i)];
            })
        .def("__setitem__", [](Step::selection& s, int i, const Atom& at){
                if (i<0){
                    i = i+static_cast<int>(s.getNat());
                }
                if ((i<0) || i>=static_cast<int>(s.getNat())){
                    throw py::index_error();
                }
                s[static_cast<size_t>(i)] = at;
        })
        .def("__len__", &Step::selection::getNat)
        .def_property_readonly("nat", &Step::selection::getNat)
        .def("__iter__", [](const Step::selection& s){return py::make_iterator(s.begin(), s.end());})
    // TYPES
        .def("getTypes", [](const Step::selection& s){
            auto oldT = s.getTypes();
            return std::vector<std::string>(oldT.begin(), oldT.end());
        })
        .def_property_readonly("ntyp", &Step::selection::getNtyp)
    // FMT
        .def("getFmt", &Step::selection::getFmt)
        .def("asFmt", py::overload_cast<AtomFmt>(&Step::selection::asFmt), "fmt"_a,
             py::return_value_policy::reference_internal)
        .def("setFmt", &Step::selection::setFmt, "fmt"_a)
    // CELL
        .def("hasCell", &Step::selection::hasCell)
        .def("enableCell", &Step::selection::enableCell, "enable"_a)
        .def("getCellDim", &Step::selection::getCellDim)
        .def("getCellVec", &Step::selection::getCellVec)
        .def("getCom", py::overload_cast<>(&Step::selection::getCom, py::const_))
        .def("getCom", py::overload_cast<AtomFmt>(&Step::selection::getCom, py::const_), "fmt"_a)
        .def("getCenter", &Step::selection::getCenter, "fmt"_a, "com"_a=false)
    // BONDS
        .def("getBonds", &Step::selection::getBonds,
             "cutfac"_a=settings.bondCutFac.val, "level"_a=settings.bondLvl.val, "update"_a=settings.bondFreq.val)
        .def("setBonds", &Step::selection::setBonds,
             "level"_a=settings.bondLvl.val, "cutfac"_a=settings.bondCutFac.val)
        .def_property_readonly("nbond", &Step::selection::getNbond)
    // SELECTION
        .def("select", [](Step::selection& s, std::string sel){return Step::selection{s.select(sel)};})
    // TODO: Modification functions
    // TODO: selection specific functions?
    ;

}
}
}
