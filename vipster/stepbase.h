#ifndef STEPBASE_H
#define STEPBASE_H

#include "atom.h"
#include "bond.h"
#include "vec.h"
#include "data.h"
#include "stepformatter.h"
#include "stepselection.h"

#include <memory>
#include <vector>
#include <set>

namespace Vipster {

    namespace detail{
        /* Helper functions to ensure step Formatters or Selections aren't nested in itself, e.g.:
         *
         * make_formatter_t<Formatter<T>> will return Formatter<T> instead of Formatter<Formatter<T>>
         *
         * elso enforce that Selection has priority over Formatter, i.e. Selection<Formatter<T>>
         */

        template<typename, template<typename...> typename>
        static constexpr bool is_instance{false};

        template<typename A, template<typename...> typename B>
        static constexpr bool is_instance<B<A>, B>{true};

        template<typename T>
        static constexpr bool is_selection = is_instance<T, Selection>;

        template<typename T>
        static constexpr bool is_formatter = is_instance<T, Formatter>;

        template<typename T>
        struct make_selection {
            using type = Selection<T>;
        };

        template<typename T>
        struct make_selection<Selection<T>> {
            using type = typename make_selection<T>::type;
        };

        template<typename T>
        using make_selection_t = typename make_selection<T>::type;

        template<typename T>
        struct make_formatter {
            using type = Formatter<T>;
        };

        template<typename T>
        struct make_formatter<Formatter<T>>{
            using type = typename make_formatter<T>::type;
        };

        template<typename T>
        struct make_formatter<Selection<T>>{
            using type = make_selection_t<typename make_formatter<T>::type>;
        };

        template<typename T>
        using make_formatter_t = typename make_formatter<T>::type;

    }

    /* Const-Base for all Step-like containers
     *
     * Implements const-interface for common state:
     * - Atoms (atomcontainer must be provided as template argument)
     * - Bonds
     * - Comment
     */
    template<typename T>
    class StepConst
    {
        template<typename U> friend class StepConst;
    public:
        virtual ~StepConst() = default;
        StepConst(const StepConst&) = default;
        StepConst(StepConst &&) = default;
        StepConst& operator=(StepConst&&) = default;
        StepConst& operator=(const StepConst&) = default;

        // Periodic table
        const PeriodicTable& getPTE() const;

        // Atoms
        using atom_source = T;
        using const_atom = typename T::const_atom;
        size_t              getNat() const noexcept;
        const_atom          operator[](size_t i) const;
        const_atom          at(size_t i) const;
        const_atom          front() const;
        const_atom          back() const;
        // TODO: remove this? breaks invariants of atomcontext
        const atom_source&  getAtoms() const noexcept;

        // Atom-iterators
        using const_iterator = typename T::const_iterator;
        using const_reverse_iterator = std::reverse_iterator<const_iterator>;
        const_iterator          begin()   const;
        const_iterator          cbegin()  const;
        const_iterator          end()     const;
        const_iterator          cend()    const;
        const_reverse_iterator  rbegin()  const;
        const_reverse_iterator  crbegin() const;
        const_reverse_iterator  rend()    const;
        const_reverse_iterator  crend()   const;

        // Selection
        using const_selection = StepConst<detail::make_selection<T>>;
        const_selection select(SelectionFilter filter) const;

        // Comment
        const std::string&  getComment() const noexcept;

        // Types
        std::set<std::string>   getTypes() const;
        size_t                  getNtyp() const;

        // Format
        using const_formatter = StepConst<detail::make_formatter_t<T>>;
        const_formatter     asFmt(AtomFmt tgt) const;
        AtomFmt             getFmt() const;

        // Topology
        const std::vector<Bond>&    getBonds() const;
        const std::vector<Overlap>& getOverlaps() const;
        std::tuple<std::vector<Angle>,
                   std::vector<Dihedral>,
                   std::vector<Dihedral>> getTopology(bool angles=true, bool dihedrals=true, bool impropers=true) const;

        // Cell
        bool        hasCell() const;
        double      getCellDim(AtomFmt fmt) const;
        const Mat&  getCellVec() const;
        Vec         getCom() const;
        Vec         getCom(AtomFmt fmt) const;
        Vec         getCenter(AtomFmt fmt, bool com=false) const;

    protected:
        StepConst(std::shared_ptr<atom_source> atoms, std::shared_ptr<BondList> bonds,
                  std::shared_ptr<std::string> comment);
        // Data
        std::shared_ptr<atom_source>    atoms;
        std::shared_ptr<BondList>       bonds;
        std::shared_ptr<std::string>    comment;
    };

    /*
     * Base for mutable Step-like containers
     *
     * Implements non-const interface and common utility functions
     *
     * Template-argument shall be atomcontainer
     */
    template<typename T>
    class StepMutable: public StepConst<T>
    {
        template<typename U> friend class StepMutable;
    public:
        virtual ~StepMutable() = default;
        StepMutable(const StepMutable &) = default;
        StepMutable(StepMutable &&) = default;
        StepMutable& operator=(const StepMutable&) = default;
        StepMutable& operator=(StepMutable&&) = default;

        // Periodic table
        PeriodicTable& getPTE();

        // Atoms
        using typename StepConst<T>::atom_source;
        using atom = typename T::atom;
        using StepConst<T>::operator[];
        atom        operator[](size_t i);
        using StepConst<T>::at;
        atom        at(size_t i);
        atom        front();
        atom        back();
        // Atom-iterators
        using iterator = typename T::iterator;
        using reverse_iterator = std::reverse_iterator<iterator>;
        using StepConst<T>::begin;
        iterator    begin();
        using StepConst<T>::end;
        iterator    end();
        using StepConst<T>::rbegin;
        reverse_iterator rbegin();
        using StepConst<T>::rend;
        reverse_iterator rend();

        // Selection
        using StepConst<T>::select;
        using selection = StepMutable<detail::make_selection_t<T>>;
        selection select(SelectionFilter filter);

        // Comment
        void setComment(const std::string& s);

        // Format
        using StepConst<T>::asFmt;
        using formatter = StepMutable<detail::make_formatter_t<T>>;
        formatter     asFmt(AtomFmt tgt);

        // Bonds
        void generateBonds(bool overlap_only=false) const;
        void addBond(size_t at1, size_t at2, DiffVec diff={}, const std::string& type="");
        void delBond(size_t idx);
        void setBondType(size_t idx, std::string type);

        // Modifier functions
        void modShift(Vec shift, double fac=1.0);
        void modRotate(double angle, Vec axis, Vec shift={0,0,0});
        void modMirror(Vec ax1, Vec ax2, Vec shift={0,0,0});

    protected:
        StepMutable(std::shared_ptr<T> atoms, std::shared_ptr<BondList> bonds,
                    std::shared_ptr<std::string> comment);

    private:
        // generate non-periodic bonds/overlaps
        template<bool overlap_only>
        void setBondsMolecule() const;

        // generate periodic bonds/overlaps
        template<bool overlap_only>
        void setBondsCell() const;
    };
}

#include "stepbase.tpp"

#endif // STEPBASE_H
