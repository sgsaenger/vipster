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

    constexpr const char* BondsAbout =
        "Bonds can be handled both manually/explicit or automatically/implicit.\n\n"
        "Per default, Molecules use automatic bond generation."
        "This can be overridden by user request, "
        "or by an input plugin if it detects explicit bonds.\n\n"
        "In manual mode the user can add/delete bonds one by one, "
        "but also request generating bonds by the same criteria as in automatic mode.\n\n"
        "The bond type will by default be defined by the types of bonded atoms. "
        "Each bond can be assigned an arbitrary string as type. "
        "This can be used to differentiate bonds when preparing inputs for forcefield calculations. "
        "Note that the strings will not be verified to be meaningful in any way."
        ;
}

#endif // BOND_H
