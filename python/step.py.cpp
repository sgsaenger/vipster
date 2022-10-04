#include "step.py.h"
#include "vipster/step.h"

#include <utility>

using namespace Vipster;

/* FIXME: Optimize assignment
 *
 * problems so far:
 * - assignment from formatter not working as intended, library problem?
 * - assignment to sub-steps not working, types mapped multiple times
 *   -> pybind/python can't deduce correct overload at RT
 * - bindings are not templated -> confined to one atom_source hierarchy
 *
 * Possible solutions:
 * - Introduce a standalone atom to allow for conversion between hierarchies
 *   (probably slow, but only when working within python)
 *   could also solve the problem of comparability
 * - bind formatter/selection within bind_step so assignments can refer to concrete instances
 *   of the whole hierarchy (better)
 *   OR
 *   remove duplicate bindings via constexpr if and detail::is_selection/is_formatter (worse)
 */
template <typename S, typename A>
void __setitem__(S& s, int i, const A& at){
    if (i<0){
        i = i+static_cast<int>(s.getNat());
    }
    if ((i<0) || i>=static_cast<int>(s.getNat())){
        throw py::index_error();
    }
    s[static_cast<size_t>(i)] = at;
}

template <typename S>
py::class_<S> bind_step(py::handle &m, std::string name){
    py::class_<StepConst<typename S::atom_source>>(m, ("__"+name+"Base").c_str());
    auto s = py::class_<S, StepConst<typename S::atom_source>>(m, name.c_str())
        .def("getPTE", &S::getPTE)
        .def_property("comment", &S::getComment, &S::setComment)
    // Atoms
        .def("__getitem__", [](S& s, int i){
                if (i<0){
                    i = i+static_cast<int>(s.getNat());
                }
                if ((i<0) || i>=static_cast<int>(s.getNat())){
                    throw py::index_error();
                }
                return s[static_cast<size_t>(i)];
            }, py::keep_alive<0, 1>())
        .def("__setitem__", __setitem__<S, typename S::atom>)
        .def("__setitem__", __setitem__<S, typename S::formatter::atom>)
        .def("__setitem__", __setitem__<S, typename S::selection::atom>)
        .def("__setitem__", __setitem__<S, typename S::formatter::selection::atom>)
        .def("__len__", &S::getNat)
        .def_property_readonly("nat", &S::getNat)
        .def("__iter__", [](S& s){return py::make_iterator(s.begin(), s.end());})
    // TYPES
        .def("getTypes", [](const S& s){
            auto oldT = s.getTypes();
            return std::vector<std::string>(oldT.begin(), oldT.end());
        })
        .def_property_readonly("ntyp", &S::getNtyp)
    // FMT
        .def_property_readonly("fmt", &S::getFmt)
        .def("asFmt", py::overload_cast<AtomFmt>(&S::asFmt), "fmt"_a)
    // CELL
        .def_property_readonly("hasCell", &S::hasCell)
        .def("getCellDim", &S::getCellDim)
        .def("getCellVec", &S::getCellVec)
        .def("getCom", py::overload_cast<>(&S::getCom, py::const_))
        .def("getCom", py::overload_cast<AtomFmt>(&S::getCom, py::const_), "fmt"_a)
        .def("getCenter", &S::getCenter, "fmt"_a, "com"_a)
    // BONDS
        .def("getBonds", [](const S& s, bool update){
                if(update){
                    s.generateBonds();
                }
                return s.getBonds();
            }, "update"_a=true)
        .def("generateBonds", &S::generateBonds, "overlap_only"_a=false)
        .def("getOverlaps", &S::getOverlaps)
        .def("getTopology", &S::getTopology, "angles"_a=true, "dihedrals"_a=true, "impropers"_a=true)
    // SELECTION
        .def("select", [](S &s, const std::string &sel){return s.select(sel);}, "selection"_a)
    // Modification functions
        .def("modShift", &S::modShift, "shift"_a, "factor"_a=1.0f)
        .def("modRotate", &S::modRotate, "angle"_a, "axis"_a, "shift"_a=Vec{})
        .def("modMirror", &S::modMirror, "axis1"_a, "axis1"_a, "shift"_a=Vec{})
    ;

    using Atom = typename S::atom;
    using _Vec = decltype(std::declval<typename S::atom>().coord);
    auto a = py::class_<Atom>(s, "Atom")
        .def_property("name", [](const Atom &a)->const std::string&{return a.name;},
                      [](Atom &a, const std::string &s){a.name = s;})
        .def_property("coord", [](const Atom &a)->const _Vec&{return a.coord;},
                      [](Atom &a, Vec c){a.coord = c;})
        .def_property("properties", [](const Atom &a)->const AtomProperties&{return a.properties;},
                      [](Atom &a, const AtomProperties &bs){a.properties = bs;})
//        .def("__eq__", [](const Atom &lhs, const Atom &rhs){return lhs == rhs;},py::is_operator())
//        .def(py::self == py::self)
//        .def(py::self != py::self)
    ;

    py::class_<_Vec>(a, "_Vec")
        .def("__getitem__", [](const _Vec& v, int i) -> Vec::value_type{
            if (i<0) {
                i += 3;
            }
            if ((i<0) || (i>=3)) {
                throw py::index_error();
            }
            return v[i];
        })
        .def("__setitem__", [](_Vec& v, int i, Vec::value_type val){
            if (i<0) {
                i += 3;
            }
            if ((i<0) || (i>=3)) {
                throw py::index_error();
            }
            v[i] = val;
        }, py::keep_alive<0, 1>())
        .def("__repr__", [name](const _Vec& v){
            return name + "::Atom::Vec["
                    + std::to_string(v[0]) + ", "
                    + std::to_string(v[1]) + ", "
                    + std::to_string(v[2]) + "]";
        })
        .def(py::self == py::self)
        .def(py::self == Vec())
    ;

    return s;
}

void Vipster::Py::Step(py::module& m){

    auto s = bind_step<Vipster::Step>(m, "Step")
        .def(py::init<AtomFmt, std::string>(), "fmt"_a=AtomFmt::Angstrom, "comment"_a="")
    // Format
        .def("setFmt", &Step::setFmt, "fmt"_a, "scale"_a=true)
    // Atoms
        .def("newAtom", [](Vipster::Step& s, std::string name, Vec coord, AtomProperties prop){
             s.newAtom(name, coord, prop);},
             "name"_a="", "coord"_a=Vec{}, "properties"_a=AtomProperties{})
        .def("newAtom", [](Vipster::Step& s, const Step::atom& at){s.newAtom(at);}, "at"_a)
        .def("newAtom", [](Vipster::Step& s, const Step::formatter::atom& at){s.newAtom(at);}, "at"_a)
        .def("newAtom", [](Vipster::Step& s, const Step::selection::atom& at){s.newAtom(at);}, "at"_a)
        .def("newAtom", [](Vipster::Step& s, const Step::selection::formatter::atom& at){s.newAtom(at);}, "at"_a)
        .def("newAtoms", [](Vipster::Step& s, size_t i){s.newAtoms(i);}, "i"_a)
        .def("newAtoms", [](Vipster::Step& s, const Vipster::Step& rhs){s.newAtoms(rhs);}, "step"_a)
        .def("newAtoms", [](Vipster::Step& s, const Vipster::Step::formatter& rhs){s.newAtoms(rhs);}, "step"_a)
        .def("newAtoms", [](Vipster::Step& s, const Vipster::Step::selection& rhs){s.newAtoms(rhs);}, "step"_a)
        .def("newAtoms", [](Vipster::Step& s, const Vipster::Step::selection::formatter& rhs){s.newAtoms(rhs);}, "step"_a)
        .def("delAtom", [](Vipster::Step& s, size_t i){s.delAtom(i);}, "i"_a)
    // CELL
        .def("enableCell", &Vipster::Step::enableCell, "val"_a)
        .def("setCellDim", &Vipster::Step::setCellDim, "cdm"_a, "fmt"_a, "scale"_a=false)
        .def("setCellVec", &Vipster::Step::setCellVec, "vec"_a, "scale"_a=false)
    // Modification functions
        .def("modWrap", &Step::modWrap)
        .def("modCrop", &Step::modCrop)
        .def("modMultiply", &Step::modMultiply, "x"_a, "y"_a, "z"_a)
        .def("modAlign", &Step::modAlign, "step_dir"_a, "target_dir"_a)
        .def("modReshape", &Step::modReshape, "newMat"_a, "newCdm"_a, "cdmFmt"_a)
    ;

    auto sf = bind_step<Vipster::Step::formatter>(s, "__Formatter");
    bind_step<Vipster::Step::selection>(s, "__Selection");
    bind_step<Vipster::Step::formatter::selection>(sf, "__Selection");
}
