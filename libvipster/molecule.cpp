#include <utility>

#include "molecule.h"

using namespace Vipster;


Molecule::Molecule(const std::string &name, unsigned long s):
    name{name},
    kpoints{}
{
    for(std::vector<StepProper>::size_type i=0;i!=s;++i){
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

std::vector<StepProper>& Molecule::newSteps(const std::vector<StepProper> &v)
{
    steps.reserve(steps.size()+v.size());
    std::vector<StepProper>::iterator pos = steps.end();
    steps.insert(steps.end(),v.begin(),v.end());
    for(;pos!=steps.end();++pos)
    {
        pos->pse = pse;
    }
    return steps;
}

StepProper& Molecule::getStep(size_t idx)
{
    return steps.at(idx);
}

const StepProper& Molecule::getStep(size_t idx) const
{
    return steps.at(idx);
}

std::vector<StepProper>& Molecule::getSteps() noexcept
{
    return steps;
}

const std::vector<StepProper>& Molecule::getSteps() const noexcept
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
