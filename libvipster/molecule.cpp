#include "molecule.h"

using namespace Vipster;


Molecule::Molecule(std::string name, ulong s, PseMap pse):
    name{name},
    pse{pse}
{
    for(std::vector<Step>::size_type i=0;i!=s;++i){
        steps.emplace_back(pse);
    }
    stepIdx = s-1;
}

void Molecule::setCellDimAll(float cdm, bool scale)
{
    for(auto s:steps){
        s.setCellDim(cdm,scale);
    }
}

Step& Molecule::curStep()
{
    return steps.at(stepIdx);
}
