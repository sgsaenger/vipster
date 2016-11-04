#ifndef MOLECULE_H
#define MOLECULE_H

#include <global.h>
#include <step.h>
#include <config.h>

namespace Vipster {
    class Molecule
    {
    public:
        Molecule(std::string name="New Molecule",ulong s=1);
//        Molecule(const Mol& m);
//        Molecule(Mol&& m);
//        Molecule(const Step&);
        std::string name;
        std::shared_ptr<PseMap> pse = std::make_shared<PseMap>(&Vipster::pse);
        void setCellDimAll(float cdm, bool scale=false, AtomFmt fmt=AtomFmt::Bohr);
        //void setCellVecAll(float v11, float v12, float v13,
        //                   float v21, float v22, float v23,
        //                   float v31, float v32, float v33,bool scale=false);
        //void setCellVecAll(t_vec v1, t_vec v2, t_vec v3,bool scale=false);
        //void setCellVecAll(std::array<t_vec,3> vec,bool scale=false);
    //private:
        std::vector<Step> steps;
        void insertStep(const Step& step);
        void insertStep(Step&& step);
    };
}
#endif // MOLECULE_H
