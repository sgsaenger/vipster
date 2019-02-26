#ifndef GUIGLOBALS_H
#define GUIGLOBALS_H

#include <type_traits>
#include <array>

namespace Vipster{
namespace Enums{
enum GuiChange{
    atoms=0x1, cell=0x2, fmt=0x4, kpoints=0x8,
    selection=0x10, definitions=0x20,
    settings=0x40,
    extra=0x80, trajec=0x100};
}
using Enums::GuiChange;
using guiChange_t = std::underlying_type<GuiChange>::type;
constexpr auto guiStepChanged = GuiChange::atoms | GuiChange::cell | GuiChange::fmt | GuiChange::selection;
constexpr auto guiMolChanged = GuiChange::kpoints;

using PBCVec = std::array<uint8_t, 3>;

}

#endif // GUIGLOBALS_H
