#ifndef ATOM_H
#define ATOM_H

#include <string>
#include <vector>
#include <bitset>

#include "vec.h"
#include "periodictable.h"

namespace Vipster{
    // supported formats
    // TODO: (dynamic?) map that stores conversion factors? unify with CdmFmt
    enum class AtomFmt { Crystal=-2, Alat, Angstrom, Bohr};
    inline bool atomFmtRelative(AtomFmt f) {return f<=AtomFmt::Alat;}
    inline bool atomFmtAbsolute(AtomFmt f) {return !atomFmtRelative(f);}
    constexpr size_t nAtFmt = 4;

    namespace Enums{
    enum AtomFlag: uint8_t {FixX, FixY, FixZ, Hidden};
    }
    using Enums::AtomFlag;
    constexpr size_t nAtFlag = 4;
    using AtomFlags = std::bitset<nAtFlag>;
    struct AtomProperties{
        double      charge;
        Vec         forces;
        AtomFlags   flags;
    };
    inline bool operator==(const AtomProperties &p1, const AtomProperties &p2){
        return std::tie(p1.charge, p1.flags, p1.forces)
               ==
               std::tie(p2.charge, p2.flags, p2.forces);
    }

    constexpr const char* AtomsAbout =
        "Atoms have the following properties:\n\n"
        "Type:\n"
        "Each atom can be assigned an arbitrary name, which will define its type. "
        "Each loaded molecule has its own periodic table where used types are saved. "
        "When a name is assigned, Vipster will try to determine a suitable type by fuzzy "
        "matching against known types in the global periodic table. "
        "This will be done by successively stripping characters from the end of the atom's name. "
        "If the name is an integer, it will match via the atom number. "
        "The properties of this type can be changed in the periodic table and will only affect "
        "the currently active molecule. "
        "Changed and/or custom types can be saved to the global table for reuse."
        "\n\n"
        "Position:\n"
        "Atomic coordinates are available in different formats:\n"
        "Absolute coordinates in either Ångström or Bohr radii, "
        "or relative coordinates scaled by the lattice constant, "
        "either relative to the cartesian axes (alat) or to the cell vectors (crystal).\n"
        "These formats can be used interchangeably. "
        "The \"active\" format will determine the base format, "
        "and will be the default for operations like script execution, "
        "if no explicit format is given."
        "\n\n"
        "Charge/Forces:\n"
        "Only used in reading and writing files."
        "\n\n"
        "Visibility:\n"
        "Flag that controls whether this atom and its bonds will be drawn."
        "\n\n"
        "Constraints:\n"
        "Flags used in some file formats that control whether this atom is allowed to move "
        "in the specified direction."
        ;
}
#endif // ATOM_H
