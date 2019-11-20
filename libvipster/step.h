#ifndef LIBVIPSTER_STEP_H
#define LIBVIPSTER_STEP_H

#include "stepmutable.h"
#include "stepsel.h"

#include <vector>
#include <memory>

namespace Vipster {

/*
 * Basic serial Atom container
 *
 * Stores atom in separate vectors
 */
struct AtomList{
    template<typename T>
    class AtomListIterator;
    using iterator = AtomListIterator<Atom>;
    using const_iterator = AtomListIterator<constAtom>;

    // Coordinates
    // one buffer per Vipster::AtomFmt
    std::array<std::vector<Vec>, nAtFmt> coordinates;
    std::array<bool, nAtFmt>             coord_changed;
    std::array<bool, nAtFmt>             coord_outdated;
    // Pointers to PeriodicTable entries
    std::vector<PeriodicTable::value_type*> elements;
    bool                                    element_changed;
    // Properties
    std::vector<AtomProperties>     properties;
    bool                            prop_changed;
    // interface
    void evaluateCache(const StepConst<AtomList>&);
    size_t getNat() const noexcept;

    template<typename T>
    class AtomListIterator: private T
    {
        template<typename U> friend class AtomListIterator;
    public:
        using difference_type = ptrdiff_t;
        using value_type = T;
        using reference = T&;
        using pointer = T*;
        using iterator_category = std::bidirectional_iterator_tag;
        AtomListIterator(const std::shared_ptr<AtomList> &atoms,
                         const std::shared_ptr<PeriodicTable> &pte,
                         AtomFmt fmt, size_t idx)
            : T{&atoms->coordinates[static_cast<uint8_t>(fmt)][idx],
                &atoms->coord_changed[static_cast<uint8_t>(fmt)],
                &atoms->elements[idx],
                &atoms->element_changed,
                &atoms->properties[idx],
                &atoms->prop_changed,
                pte.get()
            }, atoms{atoms}, fmt{fmt}, idx{idx}
        {}
        // default copy constructor
        AtomListIterator(const AtomListIterator &it) = default;
        // copy assignment just assigns the iterator, not to the atom (which has reference semantics)
        AtomListIterator& operator=(const AtomListIterator &it){
            // reset atom's pointers as it would happen in constructor
            this->coord_ptr = it.coord_ptr;
            this->coord_changed = it.coord_changed;
            this->element_ptr = it.element_ptr;
            this->element_changed = it.element_changed;
            this->prop_ptr = it.prop_ptr;
            this->prop_changed = it.prop_changed;
            this->pte = it.pte;
            // reset own properties
            atoms = it.atoms;
            fmt = it.fmt;
            idx = it.idx;
            return *this;
        }
        // allow iterator to const_iterator conversion
        template <typename U, typename R=T, typename = typename std::enable_if<std::is_same<constAtom, R>::value>::type>
        AtomListIterator(const AtomListIterator<U>& it)
            : T{*it}, atoms{it.atoms}, fmt{it.fmt}, idx{it.idx}
        {}
        AtomListIterator& operator++(){
            ++idx;
            ++(this->coord_ptr);
            ++(this->element_ptr);
            ++(this->prop_ptr);
            return *this;
        }
        AtomListIterator& operator--(){
            --idx;
            --(this->coord_ptr);
            --(this->element_ptr);
            --(this->prop_ptr);
            return *this;
        }
        AtomListIterator& operator+=(long i){
            idx += i;
            this->coord_ptr += i;
            this->element_ptr += i;
            this->prop_ptr += i;
            return *this;
        }
        AtomListIterator& operator-=(long i){
            idx -= i;
            this->coord_ptr -= i;
            this->element_ptr -= i;
            this->prop_ptr -= i;
            return *this;
        }
        AtomListIterator operator+(long i){
            AtomListIterator copy{*this};
            return copy+=i;
        }
        AtomListIterator operator-(long i){
            AtomListIterator copy{*this};
            return copy-=i;
        }
        reference operator*() const {
            // remove constness of iterator, as it is independent of constness of Atom
            return static_cast<reference>(*const_cast<AtomListIterator*>(this));
        }
        pointer operator->() const {
            return &(operator*());
        }
        bool    operator==(const AtomListIterator& rhs) const noexcept{
            return (atoms == rhs.atoms) && (fmt == rhs.fmt) && (idx == rhs.idx);
        }
        bool    operator!=(const AtomListIterator& rhs) const noexcept{
            return !(*this == rhs);
        }
        size_t getIdx() const noexcept{
            return idx;
        }
    private:
        std::shared_ptr<AtomList> atoms;
        AtomFmt fmt;
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
    Step(std::shared_ptr<PeriodicTable>, AtomFmt,
         std::shared_ptr<AtomList>, std::shared_ptr<BondList>,
         std::shared_ptr<CellData>, std::shared_ptr<std::string>);
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

    Step    asFmt(AtomFmt); // hides StepMutable::asFmt
    using StepConst<AtomList>::asFmt;

    // Atoms
    void    newAtom(std::string name="",
                    Vec coord=Vec{},
                    AtomProperties prop=AtomProperties{});
    void    newAtom(const Atom& at);
    void    newAtoms(size_t i);
    void    newAtoms(const AtomList& atoms);
    template<typename T>
    void    newAtoms(const StepConst<T>& s)
    {
        auto step = s.asFmt(at_fmt >= AtomFmt::Crystal ? AtomFmt::Bohr : at_fmt);
        const size_t nat = this->getNat() + step.getNat();
        const size_t fmt = static_cast<size_t>(at_fmt);
        // Coordinates
        AtomList& al = *this->atoms;
        al.coordinates[fmt].reserve(nat);
        al.elements.reserve(nat);
        al.properties.reserve(nat);
        al.coord_changed[fmt] = true;
        al.element_changed = true;
        al.prop_changed = true;
        if(at_fmt >= AtomFmt::Crystal){
            auto tmp = getFormatter(AtomFmt::Bohr, AtomFmt::Crystal);
            for(const auto& at: step){
                al.coordinates[fmt].push_back(tmp(at.coord));
                al.elements.push_back(&*pte->find_or_fallback(at.name));
                al.properties.push_back(at.properties);
            }
        }else{
            for(const auto& at: step){
                al.coordinates[fmt].push_back(at.coord);
                al.elements.push_back(&*pte->find_or_fallback(at.name));
                al.properties.push_back(at.properties);
            }
        }
    }
    void    delAtom(size_t i);
    template<template<typename> class T>
    void    delAtoms(SelectionBase<T, AtomList>& s)
    {
        const auto& idx = s.getIndices();
        for(auto it = idx.rbegin(); it != idx.rend(); ++it)
        {
            delAtom(it->first);
        }
        s.setFilter(SelectionFilter{});
    }

    // Cell
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
