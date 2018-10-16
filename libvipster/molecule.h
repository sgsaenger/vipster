#ifndef MOLECULE_H
#define MOLECULE_H

#include "step.h"
#include "kpoints.h"
#include <list>

namespace Vipster {

class Molecule
{
public:
    Molecule(const std::string &name="New Molecule", size_t s=1);
    template<typename T>
    Molecule(const StepConst<T>& step, std::string name="Copy of Step")
        : name{name},
          kpoints{}
    {
        steps.emplace_back(step);
        steps.back().pse = pse;
    }
    std::shared_ptr<PseMap> pse = std::make_shared<PseMap>();

    Step& newStep(const Step& step);
    Step& newStep(Step&& step={});
    Step& getStep(size_t idx);
    const Step& getStep(size_t idx) const;
    std::list<Step>& getSteps(void) noexcept;
    const std::list<Step>& getSteps(void) const noexcept;
    size_t getNstep(void) const noexcept;

    void setName(const std::string &s);
    const std::string& getName(void) const noexcept;

    void setKPoints(const KPoints &k);
    KPoints& getKPoints(void) noexcept;
    const KPoints& getKPoints(void) const noexcept;
private:
    std::list<Step> steps;
    std::string name;
    KPoints kpoints;
};
}
#endif // MOLECULE_H
