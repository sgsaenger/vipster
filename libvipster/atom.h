#ifndef ATOM_H
#define ATOM_H

#include <string>
#include "vec.h"

namespace Vipster{
    enum class AtomFmt { Bohr, Angstrom, Crystal, Alat };

    typedef std::array<bool,3> FixVec;

    struct Atom{
        std::string name = "C";
        Vec coord = {{0, 0, 0}};
        float charge = 0;
        FixVec fix = {{false, false, false}};
        bool hidden = false;
    };
    inline bool operator ==(const Atom &a1, const Atom &a2){
        return std::tie(a1.name, a1.coord, a1.charge, a1.fix, a1.hidden)
                ==
               std::tie(a2.name, a2.coord, a2.charge, a2.fix, a2.hidden);
    }
    inline bool operator !=(const Atom &a1, const Atom &a2){
        return !(a1 == a2);
    }
}

#endif // ATOM_H
