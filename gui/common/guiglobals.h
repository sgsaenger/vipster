#ifndef GUIGLOBALS_H
#define GUIGLOBALS_H

#include <type_traits>
#include <array>

namespace Vipster::GUI{

namespace Enums{
enum Change{
    // tied to active step
    atoms=0x1, // changed atomic data
    cell=0x2, // changed cell data
    fmt=0x4, // changed format
    selection=0x8, // modified selection
    definitions=0x10, // modified defined groups
    // tied to active molecule
    kpoints=0x100, // modified or swapped k-point data
    trajec=0x200, // modified or swapped non-active steps in active molecule
    // not tied to anything
    data=0x1000, // added, modified or swapped non-atomic data
    settings=0x2000, // changed global settings
    extra=0x4000, // added/removed extra data to/from visualization
};
}

using Enums::Change;
using change_t = std::underlying_type<Change>::type;
constexpr auto stepChanged = Change::atoms | Change::cell | Change::fmt | Change::selection | Change::definitions;
constexpr auto molChanged = stepChanged | Change::kpoints | Change::trajec;

using PBCVec = std::array<uint8_t, 3>;

}

#endif // GUIGLOBALS_H
