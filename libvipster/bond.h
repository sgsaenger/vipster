#ifndef BOND_H
#define BOND_H

namespace Vipster {
    struct Bond{
        std::size_t at1;
        std::size_t at2;
        float dist;
        long xdiff;
        long ydiff;
        long zdiff;
    };
}

#endif // BOND_H
