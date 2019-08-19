#ifndef BOND_H
#define BOND_H

#include <cstdint>
#include <vector>
#include <string>
#include <set>
#include <tuple>

#include "global.h"

namespace Vipster {
    enum class BondPolicy { None, Molecule, Cell };
    enum class BondFrequency { Never, Once, Always };

    struct Bond{
        std::size_t at1;
        std::size_t at2;
        float dist;
        DiffVec diff;
        const std::string* type{nullptr};
    };
    inline bool operator==(const Bond& lhs, const Bond& rhs){
        return std::tie(lhs.at1, lhs.at2, lhs.diff)
                ==
               std::tie(rhs.at1, rhs.at2, rhs.diff);
    }

    struct BondList{
        bool                    outdated{true}, setOnce{false};
        BondPolicy              level{BondPolicy::None};
        float                   cutoff_factor{-1};
        std::vector<Bond>       bonds;
        std::set<std::string>   types;
    };

}

#endif // BOND_H
