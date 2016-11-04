#include "molecule.h"

using namespace Vipster;


Molecule::Molecule(std::string name, ulong s):
    name{name}
{
    for(std::vector<Step>::size_type i=0;i!=s;++i){
        steps.emplace_back(pse);
    }
}

void Molecule::setCellDimAll(float cdm, bool scale, Fmt fmt)
{
    for(Step& s:steps){
        s.setCellDim(cdm,scale,fmt);
    }
}
