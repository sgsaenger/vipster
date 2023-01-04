#include <utility>

#include "molecule.h"

using namespace Vipster;


Molecule::Molecule(const std::string &name, size_t s):
    name{name},
    kpoints{}
{
    for(size_t i=0; i!=s; ++i){
        steps.emplace_back();
        steps.back().setPTE(pte);
    }
}

Step& Molecule::newStep(const Step &step)
{
    steps.emplace_back(step);
    steps.back().setPTE(pte);
    return steps.back();
}

Step& Molecule::newStep(Step &&step)
{
    steps.emplace_back(std::move(step));
    steps.back().setPTE(pte);
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

PeriodicTable& Molecule::getPTE()
{
    return *pte;
}

const PeriodicTable& Molecule::getPTE() const
{
    return *pte;
}

void Molecule::cleanPTE()
{
    auto in_use = this->getTypes();
    auto &pte = getPTE();
    for(auto it = pte.begin(); it != pte.end();){
        if(in_use.find(it->first) == in_use.end()){
            it = pte.erase(it);
        }else{
            ++it;
        }
    }
}

std::set<std::string> Molecule::getTypes() const
{
    std::set<std::string> tmp{};
    for(const auto& step: steps){
        for(const auto& at: step){
            tmp.insert(at.name);
        }
    }
    return tmp;
}
