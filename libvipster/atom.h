#ifndef ATOM_H
#define ATOM_H

namespace Vipster{
    struct Atom{
        std::string name;
        Vec coord;
        float charge;
        std::array<bool,3> fix;
        bool hidden;
    };
}

#endif // ATOM_H
