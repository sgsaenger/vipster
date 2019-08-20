#ifndef BOND_H
#define BOND_H

#include <cstdint>
#include <vector>
#include <string>
#include <map>
#include <tuple>

#include "global.h"

namespace Vipster {
    enum class BondMode {Manual, Automatic};

    struct Bond{
        std::size_t at1;
        std::size_t at2;
        float dist;
        DiffVec diff;
        std::pair<const std::string, ColVec>* type{nullptr};
    };
    inline bool operator==(const Bond& lhs, const Bond& rhs){
        return std::tie(lhs.at1, lhs.at2, lhs.diff)
                ==
               std::tie(rhs.at1, rhs.at2, rhs.diff);
    }

    struct BondList{
        bool                            outdated{true};
        BondMode                        mode{BondMode::Automatic};
        std::vector<Bond>               bonds;
        std::map<std::string, ColVec>   types;
    };

}

#endif // BOND_H
