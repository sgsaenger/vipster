#ifndef BOND_H
#define BOND_H

#include <cstdint>
#include <vector>
#include <string>
#include <map>
#include <tuple>

#include "global.h"

namespace Vipster {
    struct Bond{
        std::size_t at1;
        std::size_t at2;
        double      dist;
        DiffVec     diff;
        std::pair<const std::string, ColVec>* type{nullptr};
    };
    inline bool operator==(const Bond& lhs, const Bond& rhs){
        return std::tie(lhs.at1, lhs.at2, lhs.diff)
                ==
               std::tie(rhs.at1, rhs.at2, rhs.diff);
    }
    struct Overlap{
        std::size_t at1;
        std::size_t at2;
        bool        periodic; // TODO: unused
    };

    struct BondList{
        std::vector<Bond>               list;
        std::vector<Overlap>            overlaps;
        std::map<std::string, ColVec>   types;
    };

    constexpr const char* BondsAbout =
        "Bonds can be handled both manually/explicit or automatically/implicit.\n\n"
        "The default is automatic bond generation, "
        "unless the user or a file changes to explicit bonds.\n\n"
        "Automatic bonds will be determined via the involved atom type's cutoff radii "
        "and named according to their names."
        "In manual mode the user can either add, delete and modify single bonds, "
        "or explicitely request a (re-)calculation of bonds according to the automatic mechanism.\n"
        "(WARNING: This will reset all existing bonds!)\n\n"
        "Explicit bonds can be assigned an arbitrary string as type. "
        "This can be used to differentiate bonds when preparing inputs for forcefield calculations. "
        "Note that the strings will not be verified to be meaningful in any way."
        ;
}

#endif // BOND_H
