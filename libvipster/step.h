#ifndef LIBVIPSTER_STEP_H
#define LIBVIPSTER_STEP_H

#include "atom.h"
#include "cell.h"
#include "stepbase.h"
#include "vec.h"

#include <vector>
#include <memory>

namespace Vipster {

/*
 * Basic serial Atom container
 *
 * Stores atom in separate vectors
 */
struct AtomList{
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
    std::vector<float>              charges;
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
            &atoms->charges[idx],
            &atoms->properties[idx],
            &atoms->pse[idx],
            &atoms->prop_changed,
        }, atoms{atoms}, fmt{fmt}, idx{idx}
    {}
    AtomListIterator& operator++(){
        ++idx;
        ++(this->coord_ptr);
        ++(this->name_ptr);
        ++(this->charge_ptr);
        ++(this->prop_ptr);
        ++(this->pse_ptr);
        return *this;
    }
    AtomListIterator& operator+=(size_t i){
        idx += i;
        this->coord_ptr += i;
        this->name_ptr += i;
        this->charge_ptr += i;
        this->prop_ptr += i;
        this->pse_ptr += i;
        return *this;
    }
    AtomListIterator operator+(size_t i){
        AtomListIterator copy{*this};
        return copy+=i;
    }
    T&      operator*() const {
        //const-ness of iterator is separate of const-ness of atoms!
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
 * Basic Step-class
 *
 * Instantiation of Bond- and Cell-getters with AtomList as Atom-source
 */
class Step: public StepBase<Step>
{
    friend class StepSelProper;
public:
    Step& operator=(const Step& s);
    //TODO: make move-able

    // Atoms
    size_t          getNat() const noexcept;
    void            newAtom();
    void            newAtom(std::string name,
                            Vec coord=Vec{},
                            float charge=float{},
                            AtomProperties prop=AtomProperties{});
    void            newAtom(const Atom& at);
    void            newAtoms(size_t i);
    void            delAtom(long i);
    using           iterator = AtomListIterator<Atom>;
    using           constIterator = AtomListIterator<const Atom>;
    Atom            operator[](size_t i);
    iterator        begin() noexcept;
    constIterator   begin() const noexcept;
    constIterator   cbegin() const noexcept;
    iterator        end() noexcept;
    constIterator   end() const noexcept;
    constIterator   cend() const noexcept;
    // Prop-getters
    const std::vector<Vec>&             getCoords() const noexcept;
    const std::vector<std::string>&     getNames() const noexcept;
    const std::vector<PseEntry*>&       getPseEntries() const noexcept;
    const std::vector<float>&           getCharges() const noexcept;
    const std::vector<AtomProperties>&  getProperties() const noexcept;

    // Cell-setters
    void    enableCell(bool) noexcept;
    void    setCellDim(float cdm, CdmFmt at_fmt, bool scale=false);
    void    setCellVec(const Mat &vec, bool scale=false);

protected:
    Step(std::shared_ptr<PseMap>, AtomFmt, std::shared_ptr<BondList>,
         std::shared_ptr<CellData>, std::shared_ptr<std::string>,
         std::shared_ptr<AtomList>);
    Step(const Step& s);
    std::shared_ptr<AtomList>       atoms;
};

/*
 * Wrapper to access AtomFmt-cache
 */
class StepProper;
class StepFormatter: public Step
{
public:
    StepFormatter(StepProper& step, AtomFmt at_fmt);
    Step&       asFmt(AtomFmt at_fmt) override;
    const Step& asFmt(AtomFmt at_fmt) const override;
    void evaluateCache() const override;
private:
    StepProper& step;
};

//TODO: make constructible from Step
/*
 * Concrete Step
 *
 * main owner of all data
 * contains static formatter instances
 */
class StepProper: public Step
{
    friend class StepFormatter;
public:
    StepProper(std::shared_ptr<PseMap> pse = std::make_shared<PseMap>(),
               AtomFmt at_fmt=AtomFmt::Bohr,
               std::string comment="");
    StepProper(const StepProper& s);
    StepProper& operator=(const StepProper& s);
    ~StepProper() override =default;
    //TODO: make move-able

    // Format
    void            setFmt(AtomFmt at_fmt, bool scale=false);
    StepFormatter   asBohr{*this, AtomFmt::Bohr};
    StepFormatter   asAngstrom{*this, AtomFmt::Angstrom};
    StepFormatter   asCrystal{*this, AtomFmt::Crystal};
    StepFormatter   asAlat{*this, AtomFmt::Alat};
    Step&           asFmt(AtomFmt at_fmt) override;
    const Step&     asFmt(AtomFmt at_fmt) const override;

    void evaluateCache() const override;
private:
    static constexpr StepFormatter StepProper::* fmtmap[] = {
        &StepProper::asBohr,
        &StepProper::asAngstrom,
        &StepProper::asCrystal,
        &StepProper::asAlat,
    };
};
}

#endif // LIBVIPSTER_STEP_H
