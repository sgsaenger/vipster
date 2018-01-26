#ifndef STEPINTERFACE_H
#define STEPINTERFACE_H

#include "config.h"
#include "vec.h"
#include "atomref.h"
#include "bond.h"
#include "global.h"
#include <set>
#include <memory>
#include <functional>

namespace Vipster {

enum class BondLevel { None, Molecule, Cell };
enum class CdmFmt { Bohr, Angstrom };

class Step {
public:
    virtual ~Step() = default;

    // Format
    AtomFmt                 getFmt() const noexcept;

    // Atoms
    size_t                  getNat() const noexcept;
    virtual void            newAtom() = 0;
    virtual void            newAtom(const Atom &at) = 0;
    virtual void            newAtoms(size_t i) = 0;
    virtual void            delAtom(size_t idx) = 0;
    virtual AtomRef         operator[](size_t idx) = 0;
    virtual const AtomRef   operator[](size_t idx) const = 0;

    // Atom-iterator
    class iterator: private AtomRef{
    public:
        iterator(const Step* step, size_t idx);
        iterator&       operator++();
//        iterator        operator++(int);
        AtomRef&        operator*();
        const AtomRef&  operator*() const;
        bool            operator==(const iterator&) const;
        bool            operator!=(const iterator&) const;
    private:
        Step* step;
        size_t idx;
    };
    iterator        begin();
    iterator        end();
    const iterator  begin() const;
    const iterator  end() const;
    const iterator  cbegin() const;
    const iterator  cend() const;

    // Types
    std::set<std::string>   getTypes(void) const noexcept;
    size_t                  getNtyp(void) const noexcept;

    // Cell
    virtual bool            hasCell() const noexcept = 0;
    virtual void            enableCell(bool) noexcept = 0;
    virtual void            setCellDim(float cdm, CdmFmt at_fmt, bool scale=false) = 0;
    virtual float           getCellDim(CdmFmt at_fmt) const noexcept = 0;
    virtual void            setCellVec(const Mat &vec, bool scale=false) = 0;
    virtual Mat             getCellVec(void) const noexcept = 0;
    virtual Vec             getCenter(CdmFmt at_fmt, bool com=false) const noexcept = 0;

    // Comment
    virtual void                setComment(const std::string& s) = 0;
    virtual const std::string&  getComment() const noexcept = 0;

    // Bonds
    virtual const std::vector<Bond>&    getBonds(BondLevel l=BondLevel::Cell) const = 0;
    virtual const std::vector<Bond>&    getBonds(float cutfac,
                                                 BondLevel l=BondLevel::Cell) const = 0;
    virtual size_t                      getNbond(void) const noexcept = 0;

    // Public data
    std::shared_ptr<PseMap>         pse;

protected:
    Step(std::shared_ptr<PseMap> pse, AtomFmt at_fmt);
    // Atoms
    mutable bool                        at_changed{false};
    std::shared_ptr<std::vector<Vec>>   at_coord;
    virtual void                        evaluateCache() const = 0;
    // Format
    AtomFmt                 at_fmt;
    std::function<Vec(const Vec&)> getFormatter(AtomFmt source, AtomFmt target) const noexcept;
    Vec                     formatVec(Vec in, AtomFmt source, AtomFmt target) const;
    std::vector<Vec>        formatAll(std::vector<Vec> in, AtomFmt source,
                                      AtomFmt target) const;
    //TODO: benchmark:
//    std::vector<Vec>&       formatInPlace(std::vector<Vec>& in, AtomFmt source,
//                                          AtomFmt target);
};

}

#endif // STEPINTERFACE_H
