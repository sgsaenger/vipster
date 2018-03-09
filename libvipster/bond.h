#ifndef BOND_H
#define BOND_H

#include <cstdint>
#include <vector>

namespace Vipster {
    enum class BondLevel { None, Molecule, Cell };

    struct Bond{
        std::size_t at1;
        std::size_t at2;
        float dist;
        int16_t xdiff;
        int16_t ydiff;
        int16_t zdiff;
    };

    struct BondList{
        bool                outdated{true};
        BondLevel           level{BondLevel::None};
        //TODO: read from config
        float               cutoff_factor{1.1f};
        std::vector<Bond>   bonds;
    };
}

#endif // BOND_H
