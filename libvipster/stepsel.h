#ifndef LIBVIPSTER_STEPSEL_H
#define LIBVIPSTER_STEPSEL_H

#include "step.h"

namespace Vipster {

/*
 * Selection container
 *
 * contains indices of selected atoms in AtomList
 */
struct AtomSelection{
    std::vector<size_t> indices;
    std::shared_ptr<AtomList> atoms;
};

/*
 * Iterator for Atom selection
 *
 * dereferences selection-indices
 */
template<typename T>
class AtomSelIterator: private T
{
public:
    AtomSelIterator(const std::shared_ptr<AtomSelection> &selection,
                    AtomFmt fmt, size_t idx)
        : T{&selection->atoms->coordinates[static_cast<uint8_t>(fmt)][selection->indices[idx]],
            &selection->atoms->coord_changed[static_cast<uint8_t>(fmt)],
            &selection->atoms->names[selection->indices[idx]],
            &selection->atoms->name_changed,
            &selection->atoms->charges[selection->indices[idx]],
            &selection->atoms->properties[selection->indices[idx]],
            &selection->atoms->pse[selection->indices[idx]],
            &selection->atoms->prop_changed},
          selection{selection}, fmt{fmt}, idx{idx}
    {}
    AtomSelIterator& operator++(){
        ++idx;
        auto diff = selection->indices[idx] - selection->indices[idx-1];
        this->coord_ptr += diff;
        this->name_ptr += diff;
        this->charge_ptr += diff;
        this->prop_ptr += diff;
        this->pse_ptr += diff;
        return *this;
    }
    AtomSelIterator& operator+=(size_t i){
        idx += i;
        auto diff = selection->indices[idx] - selection->indices[idx-i];
        this->coord_ptr += diff;
        this->name_ptr += diff;
        this->charge_ptr += diff;
        this->prop_ptr += diff;
        this->pse_ptr += diff;
        return *this;
    }
    AtomSelIterator operator+(size_t i){
        AtomSelIterator copy{*this};
        return copy+=i;
    }
    T&  operator*() const {
        return static_cast<T&>(*const_cast<AtomSelIterator*>(this));
    }
    T*  operator->() const {
        return &(operator*());
    }
    bool    operator==(const AtomSelIterator& rhs) const noexcept{
        return (selection == rhs.selection) && (fmt == rhs.fmt) && (idx == rhs.idx);
    }
    bool    operator!=(const AtomSelIterator& rhs) const noexcept{
        return !(*this == rhs);
    }
    size_t getIdx() const noexcept{
        return idx;
    }
private:
    std::shared_ptr<AtomSelection> selection;
    AtomFmt fmt;
    size_t idx;
};

/*
 * Basic Selection-class
 *
 * Instantiation of Bond- and Cell-getters with AtomSelection as Atom-source
 */
using Criteria = std::vector<std::string>;
class StepSelection: public StepBase<StepSelection>
{
public:
    StepSelection& operator=(const StepSelection& s);
    //TODO: missing constructors?
    // TODO: Prop-getters?

    void            evaluateCache() const override;
    Criteria&       getCriteria() noexcept;
    const Criteria& getCriteria() const noexcept;

    // Atoms
    size_t          getNat() const noexcept;
    using           iterator = AtomSelIterator<Atom>;
    using           constIterator = AtomSelIterator<const Atom>;
    Atom            operator[](size_t i);
    iterator        begin() noexcept;
    constIterator   begin() const noexcept;
    constIterator   cbegin() const noexcept;
    iterator        end() noexcept;
    constIterator   end() const noexcept;
    constIterator   cend() const noexcept;
protected:
    StepSelection(std::shared_ptr<PseMap>, AtomFmt, std::shared_ptr<BondList>,
                  std::shared_ptr<CellData>, std::shared_ptr<std::string>,
                  std::shared_ptr<AtomSelection>, std::shared_ptr<Criteria>,
                  Step&);
    StepSelection(const StepSelection& s);
    std::shared_ptr<AtomSelection>  selection;
    std::shared_ptr<Criteria>       criteria;
    Step&                           step;
};

/*
 * Wrapper to access AtomFmt-cache
 */
class StepSelProper;
class StepSelFormatter: public StepSelection
{
public:
    StepSelFormatter(StepSelProper& selector, AtomFmt fmt);
    StepSelection&       asFmt(AtomFmt at_fmt) override;
    const StepSelection& asFmt(AtomFmt at_fmt) const override;
private:
    StepSelProper& selector;
};

/*
 * Concrete Selection
 *
 * has pointer to concrete Step
 * contains static formatter instances
 */
class StepSelProper: public StepSelection
{
    friend class StepSelFormatter;
public:
    StepSelProper(StepProper& step);
//    StepSelProper(StepProper& step, std::string criterion);
//    StepSelProper(StepProper& step, Criteria criteria);
    StepSelProper(const StepSelProper& s);
    StepSelProper& operator=(const StepSelProper& s);
    ~StepSelProper() override = default;

    // Format
    void                    setFmt(AtomFmt at_fmt);
    StepSelFormatter        asBohr{*this, AtomFmt::Bohr};
    StepSelFormatter        asAngstrom{*this, AtomFmt::Angstrom};
    StepSelFormatter        asCrystal{*this, AtomFmt::Crystal};
    StepSelFormatter        asAlat{*this, AtomFmt::Alat};
    StepSelection&          asFmt(AtomFmt at_fmt) override;
    const StepSelection&    asFmt(AtomFmt at_fmt) const override;
private:
    static constexpr StepSelFormatter StepSelProper::* fmtmap[] = {
        &StepSelProper::asBohr,
        &StepSelProper::asAngstrom,
        &StepSelProper::asCrystal,
        &StepSelProper::asAlat,
    };
};

}

#endif // LIBVIPSTER_STEPSEL_H
