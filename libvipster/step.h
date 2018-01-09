#ifndef STEPINTERFACE_H
#define STEPINTERFACE_H

#include "config.h"
#include "vec.h"
#include "atom.h"
#include "bond.h"
#include "global.h"
#include <set>
#include <memory>
#include <functional>

/*
 * TODO:
 *
 * need a list of pse-pointers.
 * think about noexcept?
 * Make operator[] return something different than atom
 * (should be reworked to AtomInterface or something...)
 * think about splitting bohr/angstrom from alat/crystal(/absolute?)
 */
namespace Vipster {

enum class CdmFmt{Bohr, Angstrom};

class Step {
public:
    virtual ~Step() = default;

    // Format
    AtomFmt                 getFmt() const noexcept;

    // Atoms
    size_t              getNat() const noexcept;
    virtual void        newAtom() = 0;
    virtual void        newAtom(const Atom &at) = 0;
    virtual void        newAtoms(size_t i) = 0;
    virtual void        delAtom(size_t idx) = 0;
    virtual Atom        operator[](size_t idx) = 0;
    virtual const Atom  operator[](size_t idx) const = 0;

    // Atom-iterator
    class iterator{
    public:
        iterator(const Step*, size_t);
        iterator&   operator++();
        iterator    operator++(int);
        Atom&       operator*();
        bool        operator!=(const iterator&);
    private:
        //TODO: make child of Atom-interface?
        Step*   step;
        size_t  idx;
        Atom    at;
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
    virtual void            setCellDim(float cdm, CdmFmt fmt, bool scale=false) = 0;
    virtual float           getCellDim(CdmFmt fmt) const noexcept = 0;
    virtual void            setCellVec(const Mat &vec, bool scale=false) = 0;
    virtual Mat             getCellVec(void) const noexcept = 0;
    virtual Mat             getInvVec(void) const noexcept = 0;
    virtual Vec             getCenter(CdmFmt fmt, bool com=false) const noexcept = 0;

    // Comment
    virtual void                setComment(const std::string& s) = 0;
    virtual const std::string&  getComment() const noexcept = 0;

    // Bonds
//    virtual const std::vector<Bond>&    getBonds() const = 0;
//    virtual const std::vector<Bond>&    getBonds(float cutfac) const = 0;
//    virtual const std::vector<Bond>&    getBondsCell() const = 0;
//    virtual const std::vector<Bond>&    getBondsCell(float cutfac) const = 0;
//    virtual size_t                      getNbond(void) const noexcept = 0;

    // Public data
    std::shared_ptr<PseMap>         pse;

protected:
    Step(std::shared_ptr<PseMap> pse, AtomFmt fmt);
    // Atoms
    mutable bool                        at_changed{false};
    std::shared_ptr<std::vector<Vec>>   at_coord;
    // Format
    AtomFmt                 fmt;
    std::function<Vec(Vec)> getFormatter(AtomFmt source, AtomFmt target) const noexcept;
    Vec                     formatVec(Vec in, AtomFmt source, AtomFmt target) const;
    std::vector<Vec>        formatAll(std::vector<Vec> in, AtomFmt source,
                                      AtomFmt target) const;
    std::vector<Vec>&       formatInPlace(std::vector<Vec>& in, AtomFmt source,
                                          AtomFmt target);
};

}

#endif // STEPINTERFACE_H
