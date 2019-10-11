#ifndef GUIGLOBALS_H
#define GUIGLOBALS_H

#include <type_traits>
#include <array>

namespace Vipster::GUI{
namespace Enums{
enum Change{
    atoms=0x1, cell=0x2, fmt=0x4, kpoints=0x8,
    selection=0x10, definitions=0x20,
    settings=0x40,
    extra=0x80, trajec=0x100,
    data=0x200,
};
}
using Enums::Change;
using change_t = std::underlying_type<Change>::type;
constexpr auto stepChanged = Change::atoms | Change::cell | Change::fmt | Change::selection;
constexpr auto molChanged = Change::kpoints;

using PBCVec = std::array<uint8_t, 3>;

}

#endif // GUIGLOBALS_H
