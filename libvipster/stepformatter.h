#ifndef STEPFORMATTER_H
#define STEPFORMATTER_H

#include "step.h"

namespace Vipster {

class StepProper;
class StepFormatter : public Step
{
public:
    friend class StepProper;
    StepFormatter(StepProper* s, AtomFmt fmt);
    StepFormatter(StepProper* s, const StepFormatter& rhs);
    // Atoms
    void        newAtom();
    void        newAtom(const Atom& at);
    void        newAtoms(size_t i);
    void        delAtom(size_t idx);
    Atom        operator[](size_t idx);
    const Atom  operator[](size_t idx) const;

    // Cell
    void    setCellDim(float cdm, bool scale=false);
    float   getCellDim() const noexcept;
    void    setCellVec(const Mat &vec, bool scale=false);
    Mat     getCellVec() const noexcept;
    Vec     getCenter(bool com=false) const noexcept;

    // Comment
    void                setComment(const std::string &s);
    const std::string&  getComment() const noexcept;
private:
    StepProper*   step;
    bool    at_outdated{true};
    void    evaluateChanges();
};

}

#endif // STEPFORMATTER_H
