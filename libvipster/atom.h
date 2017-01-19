#ifndef ATOM_H
#define ATOM_H

#include <vec.h>

namespace Vipster{
    enum class AtomFmt { Bohr, Angstrom, Crystal, Alat };

    typedef std::array<bool,3> FixVec;
    inline std::ostream& operator<<(std::ostream &s, const FixVec &v)
    {
        s << std::boolalpha << "FixVec["
          << v[0] << ", " << v[1] << ", " << v[2] << "]";
        return s;
    }

    struct Atom{
        std::string name;
        Vec coord;
        float charge;
        FixVec fix;
        bool hidden;
        bool operator ==(const Atom &a2) const{
            return std::tie(name, coord, charge, fix, hidden)
                    ==
                   std::tie(a2.name, a2.coord, a2.charge, a2.fix, a2.hidden);
        }
        bool operator !=(const Atom &a2) const{
            return !(*this == a2);
        }
    };
    inline std::ostream& operator<<(std::ostream &s, const Atom &a)
    {
        s << std::boolalpha
          << "Atom:\n Name: " << a.name
          << "\n Coord: [" << a.coord[0] << ", " << a.coord[1] << ", " << a.coord[2] << "]"
          << "\n Charge: " << a.charge
          << "\n Fixed: [" << a.fix[0] << ", " << a.fix[1] << ", " << a.fix[2] << "]"
          << "\n Hidden: " << a.hidden;
        return s;
    }
}

#endif // ATOM_H
