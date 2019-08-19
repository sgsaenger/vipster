#ifndef BOND_H
#define BOND_H

#include <cstdint>
#include <vector>
#include <string>
#include <set>

#include "global.h"

namespace Vipster {
    enum class BondPolicy { None, Molecule, Cell };
    enum class BondFrequency { Never, Once, Always };

    struct Bond{
        std::size_t at1;
        std::size_t at2;
        float dist;
        DiffVec diff;
        std::string* type{nullptr};
    };

    struct BondList{
        bool                    outdated{true}, setOnce{false};
        BondPolicy              level{BondPolicy::None};
        float                   cutoff_factor{-1};
        std::vector<Bond>       bonds;
        std::set<std::string>   types;
    };

}

#endif // BOND_H
