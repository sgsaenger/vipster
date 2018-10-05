#ifndef LIBVIPSTER_STEP_H
#define LIBVIPSTER_STEP_H

#include "stepmutable.h"
#include "atomcontainers.h"
#include "stepsel.h"

#include <vector>
#include <memory>

namespace Vipster {

/*
 * Main Step-class
 *
 * with AtomList as source, this is the main class to use for atom-storage
 */

class Step: public StepMutable<AtomList>
{
    template<typename T>
    friend class SelectionBase;
public:
    Step(AtomFmt at_fmt=AtomFmt::Bohr,
         std::string comment="");
    Step(const Step& s);
    template<typename T>
    Step(const StepConst<T>& s)
        : StepMutable{s.pse, s.getFmt(),
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
    ~Step() override = default;
    //TODO: make move-able

    void evaluateCache() const override;
    Step        asFmt(AtomFmt);
    using StepConst<AtomList>::asFmt;

    StepSelection   select(std::string);
    StepSelection   select(SelectionFilter);
    // TODO: move up in hierarchy?
    StepSelConst    select(std::string) const;
    StepSelConst    select(SelectionFilter) const;

    // Atoms
    void            newAtom(std::string name="",
                            Vec coord=Vec{},
                            AtomProperties prop=AtomProperties{});
    void            newAtom(const Atom& at);
    void            newAtoms(size_t i);
    template<typename T>
    void            newAtoms(const StepConst<T>& s)
    {
        const size_t oldNat = getNat();
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
    void            delAtom(size_t i);

    // Cell
    void    setCellDim(float cdm, CdmFmt at_fmt, bool scale=false);
    void    setCellVec(const Mat &vec, bool scale=false);

protected:
    Step(std::shared_ptr<PseMap>, AtomFmt,
         std::shared_ptr<AtomList>, std::shared_ptr<BondList>,
         std::shared_ptr<CellData>, std::shared_ptr<std::string>);
};

}

#endif // LIBVIPSTER_STEP_H
