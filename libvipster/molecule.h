#ifndef MOLECULE_H
#define MOLECULE_H

#include <step.h>
#include <config.h>

namespace Vipster {
    class Molecule
    {
    public:
        Molecule(std::string name="New Molecule",ulong s=1);
        std::shared_ptr<PseMap> pse = std::make_shared<PseMap>(&Vipster::pse);
        void setCellDimAll(float cdm, bool scale=false, AtomFmt fmt=AtomFmt::Bohr);
        //void setCellVecAll(float v11, float v12, float v13,
        //                   float v21, float v22, float v23,
        //                   float v31, float v32, float v33,bool scale=false);
        //void setCellVecAll(t_vec v1, t_vec v2, t_vec v3,bool scale=false);
        //void setCellVecAll(std::array<t_vec,3> vec,bool scale=false);
        void newStep(const Step& step);
        void newStep(Step&& step);
        void newSteps(const std::vector<Step> &v);
        Step& getStep(size_t idx);
        const Step& getStep(size_t idx) const;
        std::vector<Step>& getSteps(void) noexcept;
        const std::vector<Step>& getSteps(void) const noexcept;
        size_t getNstep(void) const noexcept;
        void setName(const std::string &s);
        const std::string& getName(void) const noexcept;
    private:
        std::vector<Step> steps;
        std::string name;
    };
}
#endif // MOLECULE_H
