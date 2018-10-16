#include <utility>

#include "molecule.h"

using namespace Vipster;


Molecule::Molecule(const std::string &name, size_t s):
    name{name},
    kpoints{}
{
    for(decltype(steps)::size_type i=0; i!=s; ++i){
        steps.emplace_back();
        steps.back().pse = pse;
    }
}

Step& Molecule::newStep(const Step&step)
{
    steps.push_back(step);
    steps.back().pse = pse;
    return steps.back();
}

Step& Molecule::newStep(Step&&step)
{
    steps.push_back(std::move(step));
    steps.back().pse = pse;
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
