#ifndef MOLECULE_H
#define MOLECULE_H

#include "definitions.h"
#include "step.h"
#include "config.h"

namespace Vipster {
    class Molecule
    {
    public:
        Molecule(std::string name="",ulong s=1, PseMap pse=PseMap(&Vipster::pse));
        std::string name;
        PseMap pse;
        Step& curStep();
        void setCellDimAll(float cdm, bool scale=false);
        //void setCellVecAll(float v11, float v12, float v13,
        //                   float v21, float v22, float v23,
        //                   float v31, float v32, float v33,bool scale=false);
        //void setCellVecAll(t_vec v1, t_vec v2, t_vec v3,bool scale=false);
        //void setCellVecAll(std::array<t_vec,3> vec,bool scale=false);
    //private:
        std::vector<Step> steps;
        std::vector<Step>::size_type stepIdx;
    };
}
#endif // MOLECULE_H
