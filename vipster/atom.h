#ifndef ATOM_H
#define ATOM_H

#include <string>
#include <vector>
#include <bitset>
#include <memory>
#include <functional>

#include "vec.h"
#include "periodictable.h"

namespace Vipster{

constexpr const char* AtomsAbout =
R"--(Atoms have the following properties:

Type:
Each atom can be assigned an arbitrary name, which will define its type.
Each loaded molecule has its own periodic table where used types are saved.
When a name is assigned, Vipster will try to determine a suitable type by fuzzy matching against known types in the global periodic table.
This will be done by successively stripping characters from the end of the atom's name.
If the name is an integer, it will match via the atom number.
The properties of this type can be changed in the periodic table and will only affect the currently active molecule.
Changed and/or custom types can be saved to the global table for reuse."

Position:
Atomic coordinates are available in different formats:
Absolute coordinates in either Ångström or Bohr radii, or relative coordinates scaled by the lattice constant, either relative to the cartesian axes (alat) or to the cell vectors (crystal).
These formats can be used interchangeably.
The "active" format will determine the base format, and will be the default for operations like script execution, if no explicit format is given.

Charge/Forces:
Only used in reading and writing files.

Visibility:
Flag that controls whether this atom and its bonds will be drawn.

Constraints:
Flags used in some file formats that control whether this atom is allowed to move in the specified direction.)--";

/* Atom format
 *
 * There are two relative formats for specifying coordinates:
 * - Crystal coordinates are relative to crystal cell, i.e. cellVec * cellDim
 * - Alat coordinates are relative only to the cellDim (in a regular orthogonal coordinate system)
 *
 * All other formates are absolute coordinates in a regular orthogonal coordinate system.
 * The base unit is angstrom, conversion factors shall be registered in detail::to/fromAngstrom.
 */

enum class AtomFmt {Crystal = -2, Alat, Angstrom, Bohr};
bool atomFmtRelative(AtomFmt f);
bool atomFmtAbsolute(AtomFmt f);

namespace detail{
    constexpr std::array<double, 2> toAngstrom{1, bohrrad};
    constexpr std::array<double, 2> fromAngstrom{1, invbohr};
}

/* Atom properties
 *
 * this contains all non-essential properties of an Atom,
 * i.e. everything besides position and type
 */
struct AtomProperties{
    double      charge;
    Vec         forces;

    // Bit flags
    enum Flag: uint8_t {FixX, FixY, FixZ, Hidden};
    static constexpr size_t nAtFlag = 4;
    using Flags = std::bitset<nAtFlag>;
    Flags       flags;

    bool operator==(const AtomProperties &p) const;
};

namespace detail{

    /* context for atoms
     *
     * contains the format in which atom coordinates can be given,
     * the cell-information for relative formats,
     * and the periodic table from which types will be looked up
     */
    struct AtomContext{
        struct CellData{
            bool    enabled{false};
            double  dimension{1};
            Mat     matrix{};
            Mat     inverse{};
        };

        AtomFmt fmt{AtomFmt::Angstrom};
        std::shared_ptr<PeriodicTable> pte{std::make_shared<PeriodicTable>()};
        std::shared_ptr<CellData> cell{std::make_shared<CellData>()};
    };

    // provide conversion of coordinates between contexts
    using CoordConverter = std::function<Vec(const Vec&)>;
    CoordConverter makeConverter(const AtomContext &source,
                                 const AtomContext &target);

    template<template<bool> typename AV, bool isConst>
    class AtomIterator: private AV<isConst>
    {
        // copying is templated to allow conversion to const
        template<template<bool> typename, bool> friend class AtomIterator;
    public:
        using difference_type = ptrdiff_t;
        using value_type = AV<isConst>;
        using reference = value_type&;
        using pointer = value_type*;
        using iterator_category = std::random_access_iterator_tag;

        // Constructors
        AtomIterator();
        AtomIterator(typename value_type::Source &s, size_t i);

        // Copy Constructor
        AtomIterator(const AtomIterator &it);
        template<bool t=isConst, typename = typename std::enable_if_t<t>>
        AtomIterator(const AtomIterator<AV, false> &it);

        // Copy Assignment just assigns the iterator, not to the atom (which has reference semantics)
        AtomIterator& operator=(const AtomIterator &it);
        template<bool t=isConst, typename = typename std::enable_if_t<t>>
        AtomIterator& operator=(const AtomIterator<AV, false> &it);

        // access
        reference   operator*()  const;
        pointer     operator->() const;
        reference   operator[](difference_type i);
        size_t      getIdx()     const noexcept;

        // comparison
        difference_type operator- (const AtomIterator &rhs) const;
        bool            operator==(const AtomIterator &rhs) const;
        bool            operator!=(const AtomIterator &rhs) const;
        bool            operator< (const AtomIterator &rhs) const;
        bool            operator> (const AtomIterator &rhs) const;
        bool            operator<=(const AtomIterator &rhs) const;
        bool            operator>=(const AtomIterator &rhs) const;

        // in-/decrement
        AtomIterator&   operator++();
        AtomIterator    operator++(int);
        AtomIterator&   operator+=(difference_type i);
        AtomIterator    operator+(difference_type i);
        AtomIterator&   operator--();
        AtomIterator    operator--(int);
        AtomIterator&   operator-=(difference_type i);
        AtomIterator    operator-(difference_type i);

    private:
        size_t idx;
    };
}

}

#include "atom.tpp"

#endif // ATOM_H
