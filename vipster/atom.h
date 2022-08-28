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

    namespace detail{

        /* context for atoms
         *
         * contains the format in which atom coordinates can be given,
         * the cell-information for relative formats,
         * and the periodic table from which types will be looked up
         */
        struct AtomContext{
            enum Format {Crystal = -2, Alat, Angstrom, Bohr};
            struct CellData{
                bool    enabled{false};
                double  dimension{1};
                Mat     matrix{};
                Mat     inverse{};
            };

            Format fmt{Angstrom};
            std::shared_ptr<PeriodicTable> pte{std::make_shared<PeriodicTable>()};
            std::shared_ptr<CellData> cell{std::make_shared<CellData>()};

            // conversion factors for absolute formats
            static constexpr std::array<double, 2> toAngstrom{1, bohrrad};
            static constexpr std::array<double, 2> fromAngstrom{1, invbohr};
        };

        // provide conversion of coordinates between contexts
        using CoordConverter = std::function<Vec(const Vec&)>;
        CoordConverter makeConverter(const AtomContext &source,
                                     const AtomContext &target);

        template<template<bool> typename AtomView, bool isConst>
        class AtomIterator: private AtomView<isConst>
        {
            template<template<bool> typename, bool> friend class AtomIterator;
        public:
            using difference_type = ptrdiff_t;
            using value_type = AtomView<isConst>;
            using reference = value_type&;
            using pointer = value_type*;
            using iterator_category = std::random_access_iterator_tag;

            AtomIterator()
                : value_type{}, idx{}
            {}
            AtomIterator(typename value_type::Source &s, size_t i)
                : value_type{s, i},
                  idx{i}
            {}
            // copying is templated to allow conversion to const
            // default copy constructor
            AtomIterator(const AtomIterator &it)
                : value_type{it}, idx{it.idx}
            {}
            template<bool B, bool t=isConst, typename = typename std::enable_if_t<t>>
            AtomIterator(const AtomIterator<AtomView, B> &it)
                : value_type{it}, idx{it.idx}
            {}
            // copy assignment just assigns the iterator, not to the atom (which has reference semantics)
            AtomIterator&   operator=(const AtomIterator &it){
                value_type::pointTo(it);
                idx = it.idx;
                return *this;
            }
            template<bool B, bool t=isConst, typename = typename std::enable_if_t<t>>
            AtomIterator&   operator=(const AtomIterator<AtomView, B> &it){
                value_type::pointTo(it);
                idx = it.idx;
                return *this;
            }
            // access
            reference   operator*() const {
                // remove constness of iterator, as it is independent of constness of Atom
                return static_cast<reference>(*const_cast<AtomIterator*>(this));
            }
            pointer     operator->() const {
                return &(operator*());
            }
            reference   operator[](difference_type i){
                return *operator+(i);
            }
            size_t      getIdx() const noexcept{
                return idx;
            }
            // comparison
            difference_type operator-(const AtomIterator &rhs) const{
                return getIdx() - rhs.getIdx();
            }
            bool            operator==(const AtomIterator &rhs) const{
                return (this->source == rhs.source) && (this->idx == rhs.idx);
            }
            bool            operator!=(const AtomIterator &rhs) const{
                return !operator==(rhs);
            }
            friend bool     operator< (const AtomIterator &lhs, const AtomIterator &rhs){
                return lhs.idx < rhs.idx;
            }
            friend bool     operator> (const AtomIterator &lhs, const AtomIterator &rhs){
                return lhs.idx > rhs.idx;
            }
            friend bool     operator<=(const AtomIterator &lhs, const AtomIterator &rhs){
                return lhs.idx <= rhs.idx;
            }
            friend bool     operator>=(const AtomIterator &lhs, const AtomIterator &rhs){
                return lhs.idx >= rhs.idx;
            }
            // in-/decrement
            AtomIterator&   operator++(){
                operator+=(1);
                return *this;
            }
            AtomIterator    operator++(int){
                auto copy = *this;
                operator+=(1);
                return copy;
            }
            AtomIterator&   operator+=(difference_type i){
                idx += i;
                value_type::operator+=(i);
                return *this;
            }
            AtomIterator    operator+(difference_type i){
                auto copy = *this;
                return copy+=i;
            }
            AtomIterator&   operator--(){
                operator+=(-1);
                return *this;
            }
            AtomIterator    operator--(int){
                auto copy = *this;
                operator+=(-1);
                return copy;
            }
            AtomIterator&   operator-=(difference_type i){
                operator+=(-i);
                return *this;
            }
            AtomIterator    operator-(difference_type i){
                auto copy = *this;
                return copy-=i;
            }
        private:
            size_t idx;
        };
    }

    using AtomFmt = detail::AtomContext::Format;
    bool atomFmtRelative(AtomFmt f);
    bool atomFmtAbsolute(AtomFmt f);

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
    };
    bool operator==(const AtomProperties &p1, const AtomProperties &p2);

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
