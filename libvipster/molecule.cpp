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

void Molecule::insertStep(const Step &step)
{
    steps.push_back(step);
    steps.back().pse = pse;
}

void Molecule::insertStep(Step &&step)
{
    steps.push_back(std::move(step));
    steps.back().pse = pse;
}
