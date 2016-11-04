#ifndef ATOM_H
#define ATOM_H

namespace Vipster{
    enum class AtomFmt { Bohr = 1, Angstrom = 2, Crystal = 3, Alat = 4 };

    struct Atom{
        std::string name;
        Vec coord;
        float charge;
        std::array<bool,3> fix;
        bool hidden;
    };
}

#endif // ATOM_H
