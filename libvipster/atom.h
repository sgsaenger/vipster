#ifndef ATOM_H
#define ATOM_H

#include <string>
#include <vector>
#include <type_traits>
#include "vec.h"

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

    /*
     * Access wrapper for atom properties
     */
    template<typename T, bool Name = false>
    class PropRef {
    public:
        PropRef(const T *prop, const bool *mod)
            : p_prop{const_cast<T*>(prop)},
              p_mod{const_cast<bool*>(mod)} {}
        operator T&() const { return *p_prop; }
        PropRef<T>& operator=(const T& prop){
            *p_prop = prop;
            *p_mod = true;
            return *this;
        }
        template <typename R=T, typename = std::enable_if<Name>>
        const typename R::value_type* c_str() const {return p_prop->c_str();}
        template <typename R=T, typename = std::enable_if<std::is_array<T>::value>>
        typename R::reference operator [](size_t i) {return (*p_prop)[i];}
        template <typename R=T, typename = std::enable_if<std::is_array<T>::value>>
        typename R::const_reference operator [](size_t i) const {return (*p_prop)[i];}
    private:
        T* p_prop;
        bool* p_mod;
    };
    template<typename T>
    std::ostream& operator<<(std::ostream& s, const PropRef<T>& pr)
    {
        s << (T)pr;
        return s;
    }
    template<typename T>
    std::istream& operator>>(std::istream& s, const PropRef<T>& pr)
    {
        s >> (T&)pr;
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

    bool operator ==(const Vipster::Atom &a1, const Vipster::Atom &a2);
    bool operator !=(const Vipster::Atom &a1, const Vipster::Atom &a2);
}

#endif // ATOM_H
