#include "stepformatter.h"
#include "stepproper.h"

using namespace Vipster;

StepFormatter::StepFormatter(StepProper *s, AtomFmt fmt)
    : Step{s->pse, fmt},
      step{s}
{
    at_coord = std::make_shared<std::vector<Vec>>();
}

StepFormatter::StepFormatter(StepProper *s, const StepFormatter &rhs)
    : Step{s->pse, rhs.fmt}, step{s}
{
    at_coord = std::make_shared<std::vector<Vec>>(*rhs.at_coord);
}

void StepFormatter::evaluateChanges()
{
    // let Step sort out overall situation
    step->evaluateChanges();
    // if still outdated, pull in changes from Step
    if (at_outdated) {
        //TODO: reenable
//        *at_coord = formatAll(*step->at_coord, step->fmt, fmt);
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

//TODO: below is only garbage to allow linking:

void StepFormatter::newAtom()
{
    newAtom(AtomProper{});
}

void StepFormatter::newAtom(const Atom &at)
{
    step->at_name.push_back(at.name);
    at_coord->push_back(at.coord);
    step->at_charge.push_back(at.charge);
    step->at_fix.push_back(at.fix);
    step->at_hidden.push_back(at.hidden);
    at_changed = true;
}

void StepFormatter::newAtoms(size_t i)
{
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
    step->at_name.erase(step->at_name.begin()+idx);
    at_coord->erase(at_coord->begin()+idx);
    step->at_charge.erase(step->at_charge.begin()+idx);
    step->at_fix.erase(step->at_fix.begin()+idx);
    step->at_hidden.erase(step->at_hidden.begin()+idx);
    at_changed = true;
}

Atom StepFormatter::operator[](size_t idx)
{
    return {&step->at_name[idx], &(*at_coord)[idx], &step->at_charge[idx],
            &step->at_fix[idx], &step->at_hidden[idx], &at_changed};
}

const Atom StepFormatter::operator[](size_t idx) const
{
    return {&step->at_name[idx], &(*at_coord)[idx], &step->at_charge[idx],
            &step->at_fix[idx], &step->at_hidden[idx], &at_changed};
}

void StepFormatter::setCellDim(float cdm, bool scale)
{
    //TODO
}

float StepFormatter::getCellDim() const noexcept
{
    return step->celldim;
}

void StepFormatter::setCellVec(const Mat &vec, bool scale)
{
    //TODO
}

Mat StepFormatter::getCellVec() const noexcept
{
    return step->cellvec;
}

Vec StepFormatter::getCenter(bool com) const noexcept
{
    //TODO
    return step->getCenter(com);
}
