#ifndef BOND_H
#define BOND_H

namespace Vipster {
    struct Bond{
        std::size_t at1;
        std::size_t at2;
        float dist;
        int16_t xdiff;
        int16_t ydiff;
        int16_t zdiff;
    };
}

#endif // BOND_H
