#include "stepformatter.h"
#include "stepproper.h"
#include "atomproper.h"

using namespace Vipster;

StepFormatter::StepFormatter(StepProper *s, AtomFmt fmt)
    : Step{s->pse, fmt}, step{s}
{
    at_coord = std::make_shared<std::vector<Vec>>();
}

StepFormatter::StepFormatter(StepProper *s, const StepFormatter &rhs)
    : Step{s->pse, rhs.at_fmt}, step{s}
{
    at_coord = std::make_shared<std::vector<Vec>>(*rhs.at_coord);
}

void StepFormatter::evaluateCache() const
{
    // let Step sort out overall situation
    step->evaluateCache();
    // if still outdated, pull in changes from Step
    if (at_outdated) {
        *at_coord = formatAll(*step->at_coord, step->at_fmt, at_fmt);
    }
}

void StepFormatter::setComment(const std::string &s)
{
    step->comment = s;
}

const std::string& StepFormatter::getComment() const noexcept
{
    return step->comment;
}

void StepFormatter::newAtom()
{
    newAtom(AtomProper{});
}

void StepFormatter::newAtom(const Atom &at)
{
    evaluateCache();
    step->at_name.push_back(at.name);
    at_coord->push_back(at.coord);
    step->at_charge.push_back(at.charge);
    step->at_fix.push_back(at.fix);
    step->at_hidden.push_back(at.hidden);
    at_changed = true;
}

void StepFormatter::newAtoms(size_t i)
{
    evaluateCache();
    size_t oldNat = getNat();
    size_t nat = oldNat + i;
    step->at_name.resize(nat);
    at_coord->resize(nat);
    step->at_charge.resize(nat);
    step->at_fix.resize(nat);
    step->at_hidden.resize(nat);
    at_changed = true;
}

void StepFormatter::delAtom(size_t idx)
{
    evaluateCache();
    step->at_name.erase(step->at_name.begin()+idx);
    at_coord->erase(at_coord->begin()+idx);
    step->at_charge.erase(step->at_charge.begin()+idx);
    step->at_fix.erase(step->at_fix.begin()+idx);
    step->at_hidden.erase(step->at_hidden.begin()+idx);
    at_changed = true;
}

AtomRef StepFormatter::operator[](size_t idx)
{
    evaluateCache();
    return {&step->at_name[idx], &(*at_coord)[idx], &step->at_charge[idx],
            &step->at_fix[idx], &step->at_hidden[idx], &at_changed};
}

const AtomRef StepFormatter::operator[](size_t idx) const
{
    evaluateCache();
    return {&step->at_name[idx], &(*at_coord)[idx], &step->at_charge[idx],
            &step->at_fix[idx], &step->at_hidden[idx], &at_changed};
}

void StepFormatter::setCellDim(float cdm, CdmFmt fmt, bool scale)
{
    if(!(cdm>0))throw Error("Step::setCellDim(): "
                            "cell-dimension needs to be positive");
    evaluateCache();
    if (scale && (at_fmt != AtomFmt::Crystal)) {
        float ratio = cdm / getCellDim(fmt);
        for(auto& c:*at_coord) {c *= ratio;}
    }else if (!scale && (at_fmt == AtomFmt::Crystal)) {
        float ratio = getCellDim(fmt) / cdm;
        for(auto& c:*at_coord) {c *= ratio;}
    }
    switch(fmt){
    case CdmFmt::Bohr:
        step->celldimB = cdm;
        step->celldimA = cdm*bohrrad;
        break;
    case CdmFmt::Angstrom:
        step->celldimA = cdm;
        step->celldimB = cdm*invbohr;
        break;
    }
    at_changed = true;
}

float StepFormatter::getCellDim(CdmFmt fmt) const noexcept
{
    return step->getCellDim(fmt);
}

void StepFormatter::setCellVec(const Mat &vec, bool scale)
{
    Mat inv = Mat_inv(vec);
    evaluateCache();
    if (scale) {
        if (at_fmt == AtomFmt::Crystal) {
            // do nothing but set at_changed
            step->cellvec = vec;
            step->invvec = inv;
            at_changed = true;
        } else {
            /* TODO: will cause double copy
             * as call to step->evaluateChanges
             * will set asCrystal.at_outdated unnecessarily
             * because of this.at_changed
             */
            step->asCrystal.evaluateCache();
            step->cellvec = vec;
            step->invvec = inv;
            *at_coord = formatAll(*(step->asCrystal.at_coord),
                                  AtomFmt::Crystal, at_fmt);
            at_changed = true;
        }
    } else {
        if (at_fmt != AtomFmt::Crystal) {
            // do nothing but set at_changed
            step->cellvec = vec;
            step->invvec = inv;
            at_changed = true;
        } else {
            // copy from StepProper or StepProper::asAlat
            if (at_fmt != step->at_fmt) {
                step->cellvec = vec;
                step->invvec = inv;
                *at_coord = formatAll(*(step->at_coord),
                                      step->at_fmt, AtomFmt::Crystal);
            } else {
                step->asAlat.evaluateCache();
                step->cellvec = vec;
                step->invvec = inv;
                *at_coord = formatAll(*(step->asAlat.at_coord),
                                      step->asAlat.at_fmt, AtomFmt::Crystal);
            }
            at_changed = true;
        }
    }
}

Mat StepFormatter::getCellVec() const noexcept
{
    return step->cellvec;
}

Mat StepFormatter::getInvVec() const noexcept
{
    return step->invvec;
}

Vec StepFormatter::getCenter(CdmFmt fmt, bool com) const noexcept
{
    return step->getCenter(fmt, com);
}
