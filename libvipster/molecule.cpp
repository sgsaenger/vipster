#include "molecule.h"

using namespace Vipster;


Molecule::Molecule(std::string name, ulong s):
    name{name}
{
    for(std::vector<Step>::size_type i=0;i!=s;++i){
        steps.emplace_back(pse);
    }
}

void Molecule::setCellDimAll(float cdm, bool scale, AtomFmt fmt)
{
    for(Step& s:steps){
        s.setCellDim(cdm,scale,fmt);
    }
}

void Molecule::newStep(const Step &step)
{
    steps.push_back(step);
    steps.back().pse = pse;
}

void Molecule::newStep(Step &&step)
{
    steps.push_back(std::move(step));
    steps.back().pse = pse;
}

void Molecule::newSteps(const std::vector<Step> &v)
{
    steps.reserve(steps.size()+v.size());
    steps.insert(steps.end(),v.begin(),v.end());
}

Step& Molecule::getStep(size_t idx)
{
    return steps.at(idx);
}

const Step& Molecule::getStep(size_t idx) const
{
    return steps.at(idx);
}

std::vector<Step>& Molecule::getSteps() noexcept
{
    return steps;
}

const std::vector<Step>& Molecule::getSteps() const noexcept
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
