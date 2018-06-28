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
//TODO: forces, other properties

namespace Vipster{
    namespace Enums{
    enum AtomChange{coord=0x1, name=0x2, prop=0x4, del=0x8};
    enum AtomFlag: uint8_t {FixX, FixY, FixZ, Hidden};
    }
    using Enums::AtomChange;
    using Enums::AtomFlag;
    constexpr size_t nAtFlag = 4;
    using AtomFlags = std::bitset<nAtFlag>;
    enum class AtomFmt { Bohr, Angstrom, Crystal, Alat };
    constexpr size_t nAtFmt = 4;

    struct AtomProperties{
        float charge;
        Vec forces;
        AtomFlags flags;
    };
    /*
     * Basic atom interface
     *
     * Provides access-wrappers to properties
     */
    class Atom{
    protected:
        Atom(Vec *co, bool *c_m, std::string *n, bool *n_m,
             AtomProperties *p, PseEntry** pse, bool *p_m)
            : coord_ptr{co}, coord_changed{c_m}, name_ptr{n}, name_changed{n_m},
              prop_ptr{p}, pse_ptr{pse}, prop_changed{p_m}
        {}
        Vec                     *coord_ptr;
        bool                    *coord_changed;
        std::string             *name_ptr;
        bool                    *name_changed;
        AtomProperties          *prop_ptr;
        PseEntry*               *pse_ptr;
        bool                    *prop_changed;
    public:
        PropRef<std::string, &Atom::name_ptr, &Atom::name_changed>
            name{*this};
        PropRef<Vec, &Atom::coord_ptr, &Atom::coord_changed>
            coord{*this};
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
              prop_ptr{at.prop_ptr},
              pse_ptr{at.pse_ptr},
              prop_changed{at.prop_changed} {}
        inline Atom& operator=(const Atom& at)
        {
            name = at.name;
            coord = at.coord;
            properties = at.properties;
            pse_ptr = at.pse_ptr;
            return *this;
        }
    };

    inline bool operator==(const Atom &a1, const Atom &a2){
        return std::tie(a1.name, a1.coord, a1.properties)
               ==
               std::tie(a2.name, a2.coord, a2.properties);
    }
    inline bool operator!=(const Atom &a1, const Atom &a2){
        return !(a1==a2);
    }
    inline bool operator==(const AtomProperties &p1, const AtomProperties &p2){
        return std::tie(p1.charge, p1.flags, p1.forces)
               ==
               std::tie(p2.charge, p2.flags, p2.forces);
    }

}

#endif // ATOM_H
