#ifndef LIBVIPSTER_STEP_H
#define LIBVIPSTER_STEP_H

#include "stepbase.h"

#include <vector>
#include <memory>

namespace Vipster {

/*
 * Basic serial Atom container
 *
 * Stores atom in separate vectors
 */
struct AtomList{
    template<bool isConst>
    class AtomView;
    using atom = AtomView<false>;
    using const_atom = AtomView<true>;
    template<bool isConst>
    class AtomIterator;
    using iterator = AtomIterator<false>;
    using const_iterator = AtomIterator<true>;

    AtomList(AtomFmt fmt) : fmt{fmt} {}
    AtomList(const AtomList&) = default;
    AtomList(AtomList&&) = default;
    AtomList& operator=(const AtomList&) = default;
    AtomList& operator=(AtomList&&) = default;

    // Coordinates
    std::vector<Vec> coordinates{};
    // Pointers to PeriodicTable entries
    std::vector<PeriodicTable::value_type*> elements{};
    // Properties
    std::vector<AtomProperties> properties{};
    // complete Step-interface
    AtomFmt fmt;
    size_t getNat() const noexcept {return elements.size();}

    template<bool isConst>
    class AtomView
    {
        /* Wrapper class that wraps distributed SoA data into a singular Atom-interface
         *
         * should behave like a reference
         */
    private:
        class _Vec{
        public:
            _Vec(AtomView &a): a{a} {}
            // can only be explicitely constructed to reference an Atom
            _Vec(const _Vec&) = delete;
            // assigning changes the origin
            _Vec& operator=(const Vec& rhs){
                *a.v = rhs;
                return *this;
            }
            // convert to regular Vec-reference
            operator Vec&() {return *a.v;}
            operator const Vec&() const {return *a.v;}
            // array access
            Vec::value_type& operator[](std::size_t i){
                return (*a.v)[i];
            }
            const Vec::value_type& operator[](std::size_t i) const {
                return (*a.v)[i];
            }
            // comparison
            bool operator==(const Vec& rhs) const {return *this == rhs;}
        private:
            AtomView &a;
        };
        class _Name{
        public:
            _Name(AtomView &a): a{a} {}
            // can only be explicitely constructed to reference an Atom
            _Name(const _Name&) =delete;
            // assigning a string makes this point to the
            // appropriate entry in the periodic table
            _Name& operator=(const std::string& rhs){
                *a.elem = &*a.pte->find_or_fallback(rhs);
                return *this;
            }
            // convert to regular string-reference
            operator const std::string&() const {return (*a.elem)->first;}
            // comparison
            bool operator==(const std::string& rhs) const {return (*a.elem)->first == rhs;}
            // const char* access
            const char* c_str() const{
                return (*a.elem)->first.c_str();
            }
            // I/O
            friend std::ostream& operator<<(std::ostream &s, const _Name &n){
                return s << n.c_str();
            }
            friend std::istream& operator>>(std::istream &s, _Name &n){
                s >> n;
                return s;
            }
        private:
            AtomView &a;
        };
        class _Element{
        public:
            _Element(AtomView &a): a{a} {}
            // can only be explicitely constructed to reference an Atom
            _Element(const _Element& at) =delete;
            // Allow member-access via pointer-syntax
            const Element* operator->() const{return &(*a.elem)->second;}
            // convert to reference
            operator const Element&() const {return (*a.elem)->second;}
        private:
            AtomView &a;
        };
        class _Properties{
        public:
            _Properties(AtomView &a): a{a} {}
            // can only be explicitely constructed to reference an Atom
            _Properties(const _Properties&) =delete;
            // like real reference, assigning changes the origin
            _Properties& operator=(const AtomProperties& rhs){
                *a.v = rhs;
                return *this;
            }
            // convert to regular AtomProperties-reference
            operator const AtomProperties&() const {return *a.prop;}
            // comparison
            bool operator==(const AtomProperties &rhs) const {return *a.prop == rhs;}
            // Allow pointer-syntax
            const AtomProperties* operator->() const {return a.prop;}
            AtomProperties* operator->() {return a.prop;}
        private:
            AtomView &a;
        };
   protected:
        Vec *v;
        PeriodicTable *pte;
        PeriodicTable::value_type **elem;
        AtomProperties *prop;
    public:
        // "Data", encapsulated in wrapper objects
        std::conditional_t<isConst, const _Vec, _Vec> coord{*this};
        std::conditional_t<isConst, const _Name, _Name> name{*this};
        std::conditional_t<isConst, const _Element, _Element> type{*this};
        std::conditional_t<isConst, const _Properties, _Properties> properties{*this};

        AtomView(AtomList &al,
                 PeriodicTable &pte,
                 size_t i)
            : v{&al.coordinates[i]},
              pte{&pte},
              elem{&al.elements[i]},
              prop{&al.properties[i]}
        {}
        // copying is templated to allow conversion to const
        // copy constructor creates new object pointing to same data
        AtomView(const AtomView &rhs)
            : v{rhs.v},
              pte{rhs.pte},
              elem{rhs.elem},
              prop{rhs.prop}
        {}
        template<bool B, bool t=isConst, typename = typename std::enable_if<t>::type>
        AtomView(const AtomView<B> &rhs)
            : v{rhs.v},
              pte{rhs.pte},
              elem{rhs.elem},
              prop{rhs.prop}
        {}
        // copy assignment changes data
        AtomView& operator=(const AtomView &rhs){
            coord = rhs.coord;
            name = rhs.name;
            properties = rhs.properties;
        }
        template<bool B, bool t=isConst, typename = typename std::enable_if<t>::type>
        AtomView& operator=(const AtomView<B> &rhs){
            coord = rhs.coord;
            name = rhs.name;
            properties = rhs.properties;
        }
        virtual ~AtomView() = default;
    };

    template<bool isConst>
    class AtomIterator: private AtomView<isConst>
    {
    public:
        using difference_type = ptrdiff_t;
        using value_type = AtomView<isConst>;
        using reference = value_type&;
        using pointer = value_type*;
        using iterator_category = std::random_access_iterator_tag;
        // TODO: default constructibility
        AtomIterator(AtomList &atoms,
                     PeriodicTable &pte,
                     size_t idx)
            : value_type{atoms, pte, idx},
              idx{idx}
        {}
        // copy-functions are duplicated as templates to allow for conversion to const
        // default copy constructor
        AtomIterator(const AtomIterator &it)
            : value_type{it}, idx{it.idx}
        {}
        template<bool B, bool t=isConst, typename = typename std::enable_if<t>::type>
        AtomIterator(const AtomIterator<B> &it)
            : value_type{it}, idx{it.idx}
        {}
        // copy assignment just assigns the iterator, not to the atom (which has reference semantics)
        AtomIterator&   operator=(const AtomIterator &it){
            this->v = it.v;
            this->elem = it.elem;
            this->pte = it.pte;
            this->prop = it.prop;
            idx = it.idx;
            return *this;
        }
        template<bool B, bool t=isConst, typename = typename std::enable_if<t>::type>
        AtomIterator&   operator=(const AtomIterator<B> &it){
            this->v = it.v;
            this->elem = it.elem;
            this->pte = it.pte;
            this->prop = it.prop;
            idx = it.idx;
            return *this;
        }
        // access
        reference       operator*() const {
            // remove constness of iterator, as it is independent of constness of Atom
            return static_cast<reference>(*const_cast<AtomIterator*>(this));
        }
        pointer         operator->() const {
            return &(operator*());
        }
        reference       operator[](difference_type i){
            return *operator+(i);
        }
        size_t          getIdx() const noexcept{
            return idx;
        }
        // comparison
        difference_type operator-(const AtomIterator &rhs) const noexcept{
            return idx - rhs.idx;
        }
        bool            operator==(const AtomIterator &rhs) const noexcept{
            return (this->elem == rhs.elem);
        }
        bool            operator!=(const AtomIterator &rhs) const noexcept{
            return !operator==(rhs);
        }
        friend bool     operator< (const AtomIterator &lhs, const AtomIterator &rhs){
            return lhs.idx < rhs.idx;
        }
        friend bool     operator> (const AtomIterator &lhs, const AtomIterator &rhs){
            return lhs.idx > rhs.idx;
        }
        friend bool     operator<= (const AtomIterator &lhs, const AtomIterator &rhs){
            return lhs.idx <= rhs.idx;
        }
        friend bool     operator>= (const AtomIterator &lhs, const AtomIterator &rhs){
            return lhs.idx >= rhs.idx;
        }
        // in-/decrement
        AtomIterator&   operator++(){
            operator+=(1);
            return *this;
        }
        AtomIterator&   operator--(){
            operator+=(-1);
            return *this;
        }
        AtomIterator    operator++(int){
            auto tmp = *this;
            operator+=(1);
            return tmp;
        }
        AtomIterator    operator--(int){
            auto tmp = *this;
            operator+=(-1);
            return tmp;
        }
        AtomIterator&   operator+=(difference_type i){
            idx += i;
            this->v += i;
            this->elem += i;
            this->prop += i;
            return *this;
        }
        AtomIterator&   operator-=(difference_type i){
            operator+=(-i);
            return *this;
        }
        AtomIterator    operator+(difference_type i) const {
            AtomIterator copy{*this};
            return copy+=i;
        }
        AtomIterator    operator-(difference_type i) const {
            AtomIterator copy{*this};
            return copy-=i;
        }
    private:
        size_t idx;
    };
};

/*
 * Main Step-class
 *
 * with AtomList as source, this is the main class to use for atom-storage
 */

class Step: public StepMutable<AtomList>
{
public:
    Step(AtomFmt at_fmt=AtomFmt::Bohr,
         const std::string &comment="");
    Step(const Step& s);
    Step(Step&& s);
    Step& operator=(const Step& s);
    Step& operator=(Step&& s);
    template<typename T>
    Step(const StepConst<T>& s)
        : StepMutable<AtomList>{s.pte, s.getFmt(),
                            std::make_shared<AtomList>(),
                            std::make_shared<BondList>(),
                            std::make_shared<CellData>(),
                            std::make_shared<std::string>(s.getComment())}
    {
        enableCell(s.hasCell());
        setCellDim(s.getCellDim(CdmFmt::Bohr), CdmFmt::Bohr);
        setCellVec(s.getCellVec());
        newAtoms(s);
    }

    // Atoms
    void    newAtom(std::string name="",
                    Vec coord=Vec{},
                    AtomProperties prop=AtomProperties{});
    template<typename Atom>
    void    newAtom(const Atom& at){
        AtomList &al = *atoms;
        al.coordinates.push_back(at.coord);
        al.elements.push_back(&*pte->find_or_fallback(at.name));
        al.properties.push_back(at.properties);
    }
    void    newAtoms(size_t i);
    template<typename T>
    void    newAtoms(const StepConst<T>& s)
    {
        const size_t nat = this->getNat() + s.getNat();
        const size_t fmt = static_cast<size_t>(atoms->fmt);
        // reserve space for new atoms
        AtomList& al = *this->atoms;
        al.coordinates.reserve(nat);
        al.elements.reserve(nat);
        al.properties.reserve(nat);
        if(atoms->fmt <= AtomFmt::Alat){
            // relative positioning
            if(s.getFmt() <= AtomFmt::Alat){
                // source is relative, need to convert twice
                auto step = s.asFmt(AtomFmt::Angstrom);
                auto fmt = getFormatter(AtomFmt::Angstrom, atoms->fmt);
                for(const auto& at: step){
                    al.coordinates.push_back(fmt(at.coord));
                    al.elements.push_back(&*pte->find_or_fallback(at.name));
                    al.properties.push_back(at.properties);
                }
            }else{
                // source is absolute, convert once
                auto fmt = getFormatter(s.getFmt(), atoms->fmt);
                for(const auto& at: s){
                    al.coordinates.push_back(fmt(at.coord));
                    al.elements.push_back(&*pte->find_or_fallback(at.name));
                    al.properties.push_back(at.properties);
                }
            }
        }else{
            // absolute positioning
            auto step = s.asFmt(atoms->fmt);
            for(const auto& at: step){
                al.coordinates.push_back(at.coord);
                al.elements.push_back(&*pte->find_or_fallback(at.name));
                al.properties.push_back(at.properties);
            }
        }
    }
    void    delAtom(size_t i);
//    template<template<typename> class T>
//    void    delAtoms(SelectionBase<T, AtomList>& s)
//    {
//        const auto& idx = s.getIndices();
//        for(auto it = idx.rbegin(); it != idx.rend(); ++it)
//        {
//            delAtom(it->first);
//        }
//        s.setFilter(SelectionFilter{});
//    }

    // Cell
    void enableCell(bool val) noexcept;
    void setCellDim(double cdm, CdmFmt fmt, bool scale=false);
    void setCellVec(const Mat &vec, bool scale=false);

    // Modifier functions
    void modWrap();
    void modCrop();
    void modMultiply(size_t x, size_t y, size_t z);
    void modAlign(uint8_t step_dir, uint8_t target_dir);
    void modReshape(Mat newMat, double newCdm, CdmFmt cdmFmt);
};

}
#endif // LIBVIPSTER_STEP_H
