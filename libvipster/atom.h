#ifndef ATOM_H
#define ATOM_H

#include <vec.h>

namespace Vipster{
    enum class AtomFmt { Bohr = 1, Angstrom = 2, Crystal = 3, Alat = 4 };

    struct Atom{
        std::string name;
        Vec coord;
        float charge;
        std::array<bool,3> fix;
        bool hidden;
        //let's hope this does not break the aggregate-property:
        bool operator ==(const Atom &a2) const{
            return std::tie(name, coord, charge, fix, hidden)
                    ==
                   std::tie(a2.name, a2.coord, a2.charge, a2.fix, a2.hidden);
        }
        bool operator !=(const Atom &a2) const{
            return !(*this == a2);
        }
    };
}

#endif // ATOM_H
