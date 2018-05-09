#ifndef ATOM_H
#define ATOM_H

#include <string>
#include <vector>
#include <tuple>
#include <bitset>

#include "vec.h"
#include "atompropref.h"
#include "config.h"

//TODO: track changes in properties?

namespace Vipster{
    enum class AtomFmt { Bohr, Angstrom, Crystal, Alat };
    constexpr size_t nAtFmt = 4;
    enum AtomProperty: uint8_t {FixX, FixY, FixZ, Hidden};
    constexpr size_t nAtProp = 4;
    using AtomProperties = std::bitset<nAtProp>;

    /*
     * Basic atom interface
     *
     * Provides access-wrappers to properties
     */
    class Atom{
    protected:
        Atom(Vec *co, bool *c_m, std::string *n, bool *n_m,
             float *ch, std::bitset<nAtProp> *p, PseEntry** pse, bool *p_m)
            : coord_ptr{co}, coord_changed{c_m}, name_ptr{n}, name_changed{n_m},
              charge_ptr{ch}, prop_ptr{p}, pse_ptr{pse}, prop_changed{p_m}
        {}
        Vec                     *coord_ptr;
        bool                    *coord_changed;
        std::string             *name_ptr;
        bool                    *name_changed;
        float                   *charge_ptr;
        AtomProperties          *prop_ptr;
        PseEntry*               *pse_ptr;
        bool                    *prop_changed;
    public:
        PropRef<std::string, &Atom::name_ptr, &Atom::name_changed>
            name{*this};
        PropRef<Vec, &Atom::coord_ptr, &Atom::coord_changed>
            coord{*this};
        PropRef<float, &Atom::charge_ptr, &Atom::prop_changed>
            charge{*this};
        PropRef<AtomProperties, &Atom::prop_ptr, &Atom::prop_changed>
            properties{*this};
        const PropRef<PseEntry*, &Atom::pse_ptr, nullptr>
            pse{*this};
        virtual ~Atom() = default;
        Atom() = delete;
        inline Atom(const Atom& at)
            : coord_ptr{at.coord_ptr},
              coord_changed{at.coord_changed},
              name_ptr{at.name_ptr},
              name_changed{at.name_changed},
              charge_ptr{at.charge_ptr},
              prop_ptr{at.prop_ptr},
              pse_ptr{at.pse_ptr},
              prop_changed{at.prop_changed} {}
        inline Atom& operator=(const Atom& at)
        {
            name = at.name;
            coord = at.coord;
            charge = at.charge;
            properties = at.properties;
            pse_ptr = at.pse_ptr;
            return *this;
        }
    };

    inline bool operator==(const Atom &a1, const Atom &a2){
        return std::tie(a1.name, a1.coord, a1.charge, a1.properties)
               ==
               std::tie(a2.name, a2.coord, a2.charge, a2.properties);
    }
    inline bool operator!=(const Atom &a1, const Atom &a2){
        return !(a1==a2);
    }

}

#endif // ATOM_H
