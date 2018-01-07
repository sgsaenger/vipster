#ifndef ATOM_H
#define ATOM_H

#include <string>
#include <vector>
#include <type_traits>
#include "vec.h"

/*
 * use const where possible!
 * non-const-access may trigger reevaluation of step-properties depending
 * on atom-properties (Bonds!)
 *
 * current architecture may prove problematic when only few atoms in step
 * need to be modified
 * - non-const reading triggers unnecessary mod-flag-setting
 * - but needs to be non-const for real modifications
 * workaround: const-iteration with non-const Step::operator[] access?
 * needs user-awareness -> bad!
 *
 * TODO TODO TODO
 *
 * but will stay like this for now...
 */

/*
 * TODO:
 *
 * Benchmark performance, can references reduce overhead?
 * Or virtual functions after all?
 * References prohibit copy-construction (or just assignment?), problem?
 * Move-assign/construct PropRef should be possible
 */

namespace Vipster{
    enum class AtomFmt { Bohr, Angstrom, Crystal, Alat };

    using FixVec = std::array<bool, 3>;

    template<typename T, bool Name = false>
    class PropRef {
    public:
        PropRef(const T *prop, const bool *mod)
            : p_prop{const_cast<T*>(prop)},
              p_mod{const_cast<bool*>(mod)} {}
        // TODO! EXPENSIVE! use const Atom for read-only access!
        operator T&() {*p_mod = true; return *p_prop; }
        operator const T&() const { return *p_prop; }
        PropRef<T>& operator=(const T& prop)
        {
            *p_prop = prop;
            *p_mod = true;
            return *this;
        }
        template <typename R=T, typename = std::enable_if<Name>>
        const typename R::value_type* c_str() const {return p_prop->c_str();}
        // TODO! EXPENSIVE! use const Atom for read-only access!
        template <typename R=T, typename = std::enable_if<std::is_array<T>::value>>
        typename R::reference operator [](size_t i) {*p_mod = true; return (*p_prop)[i];}
        template <typename R=T, typename = std::enable_if<std::is_array<T>::value>>
        typename R::const_reference operator [](size_t i) const {return (*p_prop)[i];}
        bool operator==(const PropRef<T>&rhs) const noexcept
        {
            return *p_prop == *(rhs.p_prop);
        }
    private:
        T* p_prop;
        bool* p_mod;
    };
    //string may fail to be converted implicitely
    inline std::ostream& operator<<(std::ostream& s, const PropRef<std::string>& pr)
    {
        s << static_cast<const std::string &>(pr);
        return s;
    }
    inline std::istream& operator>> (std::istream& s, PropRef<std::string>& pr)
    {
        s >> static_cast<std::string&>(pr);
        return s;
    }

    /*
     * Basic atom interface
     *
     * Provides access-wrappers to properties
     */
    class Atom{
    public:
        inline Atom(const std::string *n, const Vec *co, const float *ch,
                    const FixVec *f, const char *h, const bool *m)
            :name{n,m}, coord{co,m}, charge{ch,m}, fix{f,m}, hidden{h,m} {}
        virtual ~Atom() = default;
        PropRef<std::string> name;
        PropRef<Vec> coord;
        PropRef<float> charge;
        PropRef<FixVec> fix;
        PropRef<char> hidden;
    };


    /*
     * Atom that owns its data
     */
    class AtomProper: public Atom{
    public:
        inline AtomProper(std::string name="C", Vec coord={{0,0,0}}, float charge=0,
                          FixVec fix={{false,false,false}}, char hidden=false)
            : Atom{&val_name, &val_coord, &val_charge, &val_fix, &val_hidden, &mod},
              val_name{name}, val_coord{coord}, val_charge{charge}, val_fix{fix}, val_hidden{hidden} {}
        inline AtomProper(const AtomProper& rhs)
            : Atom{&val_name, &val_coord, &val_charge, &val_fix, &val_hidden, &mod},
              val_name{rhs.name}, val_coord(rhs.coord), val_charge{rhs.charge},
              val_fix(rhs.fix), val_hidden{rhs.hidden} {}
        inline AtomProper& operator=(const AtomProper& rhs){
            val_name = rhs.val_name;
            val_coord = rhs.val_coord;
            val_charge = rhs.val_charge;
            val_fix = rhs.val_fix;
            val_hidden = rhs.val_hidden;
            return *this;
        }
    private:
        std::string val_name;
        Vec val_coord;
        float val_charge;
        FixVec val_fix;
        char val_hidden;
        bool mod;
    };

    //TODO: maybe ignore fix/hidden?
    inline bool operator==(const Atom &a1, const Atom &a2) {
        return std::tie(a1.name, a1.coord, a1.charge, a1.fix, a1.hidden)
               ==
               std::tie(a2.name, a2.coord, a2.charge, a2.fix, a2.hidden);
    }
    inline bool operator!=(const Atom &a1, const Atom &a2) {
        return !(a1==a2);
    }
}

#endif // ATOM_H
