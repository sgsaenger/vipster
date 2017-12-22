#ifndef STEP_H
#define STEP_H

#include "config.h"
#include "vec.h"
#include "atom.h"
#include "bond.h"
#include "global.h"
#include <set>
#include <memory>

/*
 * TODO:
 *
 * need a list of pse-pointers.
 *
 */

namespace Vipster{
class Step{
public:
    virtual ~Step() = default;
    // Atoms
    size_t          getNat() const noexcept;
    void            newAtom(const Atom &at=AtomProper{});
    void            newAtoms(size_t i);
    void            delAtom(size_t idx);
    Atom            operator[](size_t idx);
    const Atom      operator[](size_t idx) const;

    // Atom-iterator
    class iterator{
    public:
        iterator(const Step*, size_t);
        iterator&   operator++();
        iterator    operator++(int);
        Atom&       operator*();
        bool        operator!=(const iterator&);
    private:
        //TODO: weak_ptr to properties?
        Step*   step;
        size_t  idx;
        Atom    at;
    };
    iterator                begin();
    const iterator          begin() const;
    iterator                end();
    const iterator          end() const;

    // Types
    std::set<std::string>   getTypes(void)const noexcept;
    size_t                  getNtyp(void) const noexcept;

    // Format
    AtomFmt                 getFmt() const noexcept;

    // Cell
    virtual void            setCellDim(float cdm, bool scale=false) = 0;
    virtual float           getCellDim() const noexcept = 0;
            void            setCellVec(const Mat &vec, bool scale=false);
            Mat             getCellVec(void) const noexcept;
            Vec             getCenter(bool com=false) const noexcept;

    // Bonds
    const std::vector<Bond>&    getBonds() const;
    const std::vector<Bond>&    getBonds(float cutfac) const;
    const std::vector<Bond>&    getBondsCell() const;
    const std::vector<Bond>&    getBondsCell(float cutfac) const;
    size_t                      getNbond(void) const noexcept;

    // Public data
    std::shared_ptr<PseMap>         pse;
    std::shared_ptr<std::string>    comment;

protected:
    enum class BondLevel { None, Molecule, Cell };
    Step(const std::shared_ptr<PseMap>& pse,
         const std::shared_ptr<std::string>& comment,
         const std::shared_ptr<std::vector<std::string>>& at_name,
         const std::shared_ptr<std::vector<Vec>>& at_coord,
         const std::shared_ptr<std::vector<float>>& at_charge,
         const std::shared_ptr<std::vector<FixVec>>& at_fix,
         const std::shared_ptr<std::vector<char>>& at_hidden,
         const std::shared_ptr<float>& celldim,
         const std::shared_ptr<Mat>& cellvec,
         const std::shared_ptr<Mat>& invvec,
         const std::shared_ptr<BondLevel>& bonds_level,
         const std::shared_ptr<float>& bondcut_factor,
         const std::shared_ptr<std::vector<Bond>>& bonds);
    // Atoms
    bool                                        at_changed{false};
    std::shared_ptr<std::vector<std::string>>   at_name;
    std::shared_ptr<std::vector<Vec>>           at_coord;
    std::shared_ptr<std::vector<float>>         at_charge;
    std::shared_ptr<std::vector<FixVec>>        at_fix;
    std::shared_ptr<std::vector<char>>          at_hidden;
    // Format
    AtomFmt             fmt;
    Vec                 formatVec(Vec in, AtomFmt source, AtomFmt target) const;
    std::vector<Vec>    formatAll(const std::vector<Vec>& in, AtomFmt source,
                                  AtomFmt target) const;
    std::vector<Vec>&   formatInPlace(std::vector<Vec>& in, AtomFmt source,
                                      AtomFmt target);
    // Cell
    std::shared_ptr<float>  celldim;
    std::shared_ptr<Mat>    cellvec;
    std::shared_ptr<Mat>    invvec;
    // Bonds
    bool                                bonds_outdated{false};
    std::shared_ptr<BondLevel>          bonds_level;
    std::shared_ptr<float>              bondcut_factor;
    std::shared_ptr<std::vector<Bond>>  bonds;
    void                                setBonds(float cutfac) const;
    void                                setBondsCell(float cutfac) const;
    void                                checkBond(std::size_t i, std::size_t j,
                                                  float cutfac, Vec dist,
                                                  std::array<int, 3> offset) const;
};

class StepProper: public Step
{
public:
    StepProper();
    StepProper(const std::shared_ptr<PseMap> &pse);
//    Step(const Step&);
//    Step& operator=(const Step&);

    // FMT
    void                        setFmt(AtomFmt fmt, bool scale=false);
//    class StepFormatter: public Step{
//    public:
//        StepFormatter(StepProper *, AtomFmt);
////        Atom operator[](size_t idx);
//    private:
//        StepProper    *step;
//    };
//    const StepFormatter asAlat{this, AtomFmt::Alat};
//    const StepFormatter asAngstrom{this, AtomFmt::Angstrom};
//    const StepFormatter asBohr{this, AtomFmt::Bohr};
//    const StepFormatter asCrystal{this, AtomFmt::Crystal};
    // CELL
    void                        setCellDim(float cdm, bool scale=false);
    float                       getCellDim() const noexcept;
};

}
#endif // STEP_H
