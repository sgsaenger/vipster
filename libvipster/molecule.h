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
        pte = std::make_shared<PeriodicTable>(step.getPTE());
        steps.emplace_back();
        steps.back().atoms->ctxt.pte = pte;
        steps.back() = step;
    }

    std::string name;
    KPoints kpoints;
    // TODO: make sure that this has the appropriate root!
    std::shared_ptr<PeriodicTable> pte = std::make_shared<PeriodicTable>();

    Step& newStep(const Step& step);
    Step& newStep(Step&& step={});
    Step& getStep(size_t idx);
    const Step& getStep(size_t idx) const;
    std::list<Step>& getSteps(void) noexcept;
    const std::list<Step>& getSteps(void) const noexcept;
    size_t getNstep(void) const noexcept;

private:
    std::list<Step> steps;
};
}
#endif // MOLECULE_H
