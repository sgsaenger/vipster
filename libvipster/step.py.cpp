#include "pyvipster.h"
#include "step.h"

// TODO: think about const-correctness!
// can wait until "view"-interface to running libvipster instance is implemented

namespace Vipster::Py{
template <typename S>
py::class_<S> bind_step(py::handle &m, std::string name){
    py::class_<StepConst<typename S::source>>(m, ("__"+name+"Base__").c_str());
    auto cl = py::class_<S, StepConst<typename S::source>>(m, name.c_str())
        .def_readonly("pse", &S::pse)
        .def_property("comment", &S::getComment, &S::setComment)
    // Atoms
        .def("__getitem__", [](S& s, int i){
                s.evaluateCache();
                if (i<0){
                    i = i+static_cast<int>(s.getNat());
                }
                if ((i<0) || i>=static_cast<int>(s.getNat())){
                    throw py::index_error();
                }
                return s[static_cast<size_t>(i)];
            }, py::keep_alive<0, 1>())
        .def("__setitem__", [](S& s, int i, const Atom& at){
                s.evaluateCache();
                if (i<0){
                    i = i+static_cast<int>(s.getNat());
                }
                if ((i<0) || i>=static_cast<int>(s.getNat())){
                    throw py::index_error();
                }
                s[static_cast<size_t>(i)] = at;
        })
        .def("__len__", [](const S& s){s.evaluateCache(); return s.getNat();})
        .def_property_readonly("nat", [](const S& s){s.evaluateCache(); return s.getNat();})
        .def("__iter__", [](S& s){s.evaluateCache(); return py::make_iterator(s.begin(), s.end());})
    // TYPES
        .def("getTypes", [](const S& s){
            s.evaluateCache();
            auto oldT = s.getTypes();
            return std::vector<std::string>(oldT.begin(), oldT.end());
        })
        .def_property_readonly("ntyp", [](const S& s){s.evaluateCache(); return s.getNtyp();})
    // FMT
        .def_property("fmt", &S::getFmt, &S::setFmt)
        .def("asFmt", [](S& s, AtomFmt fmt){s.evaluateCache(); return s.asFmt(fmt);}, "fmt"_a)
    // CELL
        .def_property("hasCell", &S::hasCell, &S::enableCell)
        .def("getCellDim", &S::getCellDim)
        .def("getCellVec", &S::getCellVec)
        .def("getCom", [](const S& s){s.evaluateCache(); return s.getCom();})
        .def("getCom", [](const S& s, AtomFmt fmt){s.evaluateCache(); return s.getCom(fmt);}, "fmt"_a)
        .def("getCenter", [](const S& s, CdmFmt fmt, bool com){s.evaluateCache(); return s.getCenter(fmt, com);},
             "fmt"_a, "com"_a=false)
    // BONDS
        .def("getBonds", [](const S& s, float c, BondLevel l, BondFrequency u){
                s.evaluateCache();
                return s.getBonds(c,l,u);
            }, "cutfac"_a=settings.bondCutFac.val, "level"_a=settings.bondLvl.val, "update"_a=settings.bondFreq.val)
        .def("setBonds", &S::setBonds,
             "level"_a=settings.bondLvl.val, "cutfac"_a=settings.bondCutFac.val)
        .def_property_readonly("nbond", [](const S& s){s.evaluateCache(); return s.getNbond();})
    // SELECTION
        .def("select", [](S& s, std::string sel)->Step::selection{s.evaluateCache(); return s.select(sel);}, "sel"_a)
    ;
    return std::move(cl);
}

void Step(py::module& m){
    py::enum_<CdmFmt>(m, "CellDimFmt")
        .value("Bohr", CdmFmt::Bohr)
        .value("Angstrom", CdmFmt::Angstrom)
    ;

    auto s = bind_step<Vipster::Step>(m, "Step")
        .def(py::init<AtomFmt, std::string>(), "fmt"_a=AtomFmt::Bohr, "comment"_a="")
    // Atoms
        .def("newAtom", [](Vipster::Step& s, std::string name, Vec coord, AtomProperties prop){
             s.evaluateCache();
             s.newAtom(name, coord, prop);},
             "name"_a="", "coord"_a=Vec{}, "properties"_a=AtomProperties{})
        .def("newAtom", [](Vipster::Step& s, const Atom& at){ s.evaluateCache(); s.newAtom(at);}, "at"_a)
        .def("delAtom", [](Vipster::Step& s, size_t i){s.evaluateCache(); s.delAtom(i);}, "i"_a)
    // CELL
        .def("setCellDim", &Vipster::Step::setCellDim, "cdm"_a, "fmt"_a, "scale"_a=false)
        .def("setCellVec", &Vipster::Step::setCellVec, "vec"_a, "scale"_a=false)
    // TODO: Modification functions
    ;

    bind_step<Vipster::Step::selection>(s, "Selection");
    // TODO: Modification functions
    // TODO: selection specific functions?

}
}
