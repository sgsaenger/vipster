#include <utility>

#include "molecule.h"

using namespace Vipster;


Molecule::Molecule(const std::string &name, size_t s):
    name{name},
    kpoints{}
{
    for(decltype(steps)::size_type i=0; i!=s; ++i){
        steps.emplace_back();
        steps.back().atoms->ctxt.pte = pte;
    }
}

Step& Molecule::newStep(const Step &step)
{
    steps.push_back(step);
    steps.back().atoms->ctxt.pte = pte;
    return steps.back();
}

Step& Molecule::newStep(Step &&step)
{
    steps.push_back(std::move(step));
    steps.back().atoms->ctxt.pte = pte;
    return steps.back();
}

Step& Molecule::getStep(size_t idx)
{
    return *std::next(steps.begin(), idx);
}

const Step& Molecule::getStep(size_t idx) const
{
    return *std::next(steps.begin(), idx);
}

std::list<Step>& Molecule::getSteps() noexcept
{
    return steps;
}

const std::list<Step>& Molecule::getSteps() const noexcept
{
    return steps;
}

size_t Molecule::getNstep() const noexcept
{
    return steps.size();
}
