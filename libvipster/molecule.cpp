#include <utility>

#include "molecule.h"

using namespace Vipster;


Molecule::Molecule(const std::string &name, unsigned long s):
    name{name},
    kpoints{}
{
    for(decltype(steps)::size_type i=0; i!=s; ++i){
        steps.emplace_back(pse);
    }
}

StepProper& Molecule::newStep(const StepProper &step)
{
    steps.push_back(step);
    steps.back().pse = pse;
    return steps.back();
}

StepProper& Molecule::newStep(StepProper &&step)
{
    steps.push_back(std::move(step));
    steps.back().pse = pse;
    return steps.back();
}

StepProper& Molecule::getStep(size_t idx)
{
    return *std::next(steps.begin(), idx);
}

const StepProper& Molecule::getStep(size_t idx) const
{
    return *std::next(steps.begin(), idx);
}

std::list<StepProper>& Molecule::getSteps() noexcept
{
    return steps;
}

const std::list<StepProper>& Molecule::getSteps() const noexcept
{
    return steps;
}

size_t Molecule::getNstep() const noexcept
{
    return steps.size();
}

void Molecule::setName(const std::string &s)
{
    name = s;
}

const std::string& Molecule::getName(void)const noexcept
{
    return name;
}

void Molecule::setKPoints(const KPoints &k)
{
    kpoints = k;
}

const KPoints& Molecule::getKPoints() const noexcept
{
    return kpoints;
}

KPoints& Molecule::getKPoints() noexcept
{
    return kpoints;
}
