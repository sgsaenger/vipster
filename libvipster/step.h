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
template<typename T>
class AtomListIterator;
struct AtomList{
    using iterator = AtomListIterator<Atom>;
    using constIterator = AtomListIterator<const Atom>;
    // Coordinates
    // one buffer per Vipster::AtomFmt
    std::array<std::vector<Vec>, nAtFmt> coordinates;
    std::array<bool, nAtFmt>             coord_changed;
    std::array<bool, nAtFmt>             coord_outdated;
    // Names (synced with type-pointers)
    std::vector<std::string>        names;
    bool                            name_changed;
    std::vector<PseEntry*>          pse;
    // Properties
    std::vector<AtomProperties>     properties;
    bool                            prop_changed;
};

/*
 * Iterator for serial Atom container
 */
template<typename T>
class AtomListIterator: private T
{
public:
    AtomListIterator(const std::shared_ptr<AtomList> &atoms,
                     AtomFmt fmt, size_t idx)
        : T{&atoms->coordinates[static_cast<uint8_t>(fmt)][idx],
            &atoms->coord_changed[static_cast<uint8_t>(fmt)],
            &atoms->names[idx],
            &atoms->name_changed,
            &atoms->properties[idx],
            &atoms->pse[idx],
            &atoms->prop_changed,
        }, atoms{atoms}, fmt{fmt}, idx{idx}
    {}
    AtomListIterator& operator++(){
        ++idx;
        ++(this->coord_ptr);
        ++(this->name_ptr);
        ++(this->prop_ptr);
        ++(this->pse_ptr);
        return *this;
    }
    AtomListIterator& operator+=(size_t i){
        idx += i;
        this->coord_ptr += i;
        this->name_ptr += i;
        this->prop_ptr += i;
        this->pse_ptr += i;
        return *this;
    }
    AtomListIterator operator+(size_t i){
        AtomListIterator copy{*this};
        return copy+=i;
    }
    T&      operator*() const {
        // remove constness of iterator, as it is independent of constness of Atom
        return static_cast<T&>(*const_cast<AtomListIterator*>(this));
    }
    T*      operator->() const {
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

/*
 * Main Step-class
 *
 * with AtomList as source, this is the main class to use for atom-storage
 */

class Step: public StepMutable<AtomList>
{
public:
    Step(AtomFmt at_fmt=AtomFmt::Bohr,
         std::string comment="");
    Step(std::shared_ptr<PseMap>, AtomFmt,
         std::shared_ptr<AtomList>, std::shared_ptr<BondList>,
         std::shared_ptr<CellData>, std::shared_ptr<std::string>);
    Step(const Step& s);
    template<typename T>
    Step(const StepConst<T>& s)
        : StepMutable<AtomList>{s.pse, s.getFmt(),
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
//    Step(Step&& s);
    Step& operator=(const Step& s);
//    Step& operator=(Step&& s);

    Step    asFmt(AtomFmt); // hides StepMutable::asFmt
    using StepConst<AtomList>::asFmt;

    StepSelection   select(std::string);
    StepSelection   select(SelectionFilter);
    StepSelConst    select(std::string) const;
    StepSelConst    select(SelectionFilter) const;

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
        const size_t oldNat = this->getNat();
        const size_t nat = oldNat + s.getNat();
        // Coordinates
        AtomList& al = *this->atoms;
        al.coordinates[static_cast<size_t>(at_fmt)].reserve(nat);
        al.names.reserve(nat);
        al.properties.reserve(nat);
        al.pse.reserve(nat);
        al.coord_changed[static_cast<size_t>(at_fmt)] = true;
        al.name_changed = true;
        al.prop_changed = true;
        auto i = oldNat;
        for(const auto& at: s){
            al.coordinates[static_cast<size_t>(at_fmt)][i] = at.coord;
            al.names[i] = at.name;
            al.pse[i] = &(*pse)[at.name];
            al.properties[i] = at.properties;
        }
    }
    void    delAtom(size_t i);

    // Cell
    void    setCellDim(float cdm, CdmFmt at_fmt, bool scale=false);
    void    setCellVec(const Mat &vec, bool scale=false);

    // Modifier functions
    void modWrap();
    void modCrop();
    void modMultiply(size_t x, size_t y, size_t z);
    void modAlign(uint8_t step_dir, uint8_t target_dir);
    void modReshape(Mat newMat, float newCdm, CdmFmt cdmFmt);
};

template<>
size_t StepConst<AtomList>::getNat() const noexcept;

template<>
void StepConst<AtomList>::evaluateCache() const;

}

#endif // LIBVIPSTER_STEP_H
