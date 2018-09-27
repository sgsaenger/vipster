#ifndef LIBVIPSTER_STEP_H
#define LIBVIPSTER_STEP_H

#include "stepbase.h"
#include "atomcontainers.h"
#include "stepsel.h"

#include <vector>
#include <memory>

namespace Vipster {

/*
 * Basic Step-class
 *
 * Instantiation of Bond- and Cell-getters with AtomList as Atom-source
 */

class Step: public StepBase<Step>
{
    template<typename T>
    friend class SelectionBase;
public:
    Step& operator=(const Step& s);
    ~Step() override = default;
    //TODO: make move-able

    StepSelection   select(std::string);
    StepSelection   select(SelectionFilter);
    StepSelConst    select(std::string) const;
    StepSelConst    select(SelectionFilter) const;

    // Atoms
    size_t          getNat() const noexcept;
    void            newAtom();
    void            newAtom(std::string name,
                            Vec coord=Vec{},
                            AtomProperties prop=AtomProperties{});
    void            newAtom(const Atom& at);
    void            newAtoms(size_t i);
    void            newAtoms(const AtomList& atoms);
    void            delAtom(size_t i);
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
    const AtomList&                     getAtoms() const noexcept;
    const std::vector<Vec>&             getCoords() const noexcept;
    const std::vector<std::string>&     getNames() const noexcept;
    const std::vector<PseEntry*>&       getPseEntries() const noexcept;
    const std::vector<AtomProperties>&  getProperties() const noexcept;

    // Cell-setters
    void    enableCell(bool) noexcept;
    void    setCellDim(float cdm, CdmFmt at_fmt, bool scale=false);
    void    setCellVec(const Mat &vec, bool scale=false);
    Vec     getCenter(CdmFmt fmt, bool com=false) const noexcept override;

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
    ~StepFormatter() override = default;
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
    ~StepProper() override = default;
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
