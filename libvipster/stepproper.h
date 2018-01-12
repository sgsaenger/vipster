#ifndef STEP_H
#define STEP_H

#include "step.h"
#include "stepformatter.h"

namespace Vipster{

class StepFormatter;
class StepProper: public Step
{
public:
    friend class StepFormatter;
    //TODO: reorder arguments
    StepProper(std::shared_ptr<PseMap> pse = std::make_shared<PseMap>(),
         AtomFmt at_fmt=AtomFmt::Bohr,
         std::string comment="");
    StepProper(const StepProper& s);
    StepProper& operator=(const StepProper& s);

    // Format
    void            setFmt(AtomFmt at_fmt, bool scale=false);
    StepFormatter   asAlat;
    StepFormatter   asAngstrom;
    StepFormatter   asBohr;
    StepFormatter   asCrystal;
    StepFormatter&  asFmt(AtomFmt at_fmt);
    static constexpr StepFormatter StepProper::* fmtmap[] = {
        &StepProper::asBohr,
        &StepProper::asAngstrom,
        &StepProper::asCrystal,
        &StepProper::asAlat,
    };

    // Atoms
    void            newAtom();
    void            newAtom(const Atom& at);
    void            newAtoms(size_t i);
    void            delAtom(size_t idx);
    AtomRef         operator[](size_t idx);
    const AtomRef   operator[](size_t idx) const;

    // Cell
    void    setCellDim(float cdm, CdmFmt at_fmt, bool scale=false);
    float   getCellDim(CdmFmt at_fmt) const noexcept;
    void    setCellVec(const Mat &vec, bool scale=false);
    Mat     getCellVec() const noexcept;
    Mat     getInvVec(void) const noexcept;
    Vec     getCenter(CdmFmt at_fmt, bool com=false) const noexcept;

    // Comment
    void                setComment(const std::string &s);
    const std::string&  getComment() const noexcept;
protected:
    void                        evaluateCache() const;
private:
    //Atoms
    std::vector<std::string>    at_name;
    std::vector<float>          at_charge;
    std::vector<FixVec>         at_fix;
    std::vector<char>           at_hidden;
    //Cell
    bool        cell_enabled{true};
    float       celldimB{1};
    float       celldimA{bohrrad};
    Mat         cellvec{{{{1,0,0}}, {{0,1,0}}, {{0,0,1}}}};
    Mat         invvec{{{{1,0,0}}, {{0,1,0}}, {{0,0,1}}}};
    //Comment
    std::string comment;
    // Bonds
    enum class BondLevel { None, Molecule, Cell };
    mutable bool                bonds_outdated{false};
    mutable BondLevel           bonds_level;
    mutable float               bondcut_factor;
    mutable std::vector<Bond>   bonds;
    void                setBonds(BondLevel l, float cutfac) const;
    void                setBondsMol(float cutfac) const;
    void                setBondsCell(float cutfac) const;
    void                checkBond(std::size_t i, std::size_t j,
                                  float cutfac, const Vec& dist,
                                  const std::array<int,3> &offset) const;
};

}
#endif // STEP_H
