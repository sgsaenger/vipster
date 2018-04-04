#ifndef MOLECULE_H
#define MOLECULE_H

#include "step.h"
#include "config.h"
#include "kpoints.h"

namespace Vipster {

class Molecule
{
public:
    Molecule(const std::string &name="New Molecule",unsigned long s=1);
    std::shared_ptr<PseMap> pse = std::make_shared<PseMap>();


    StepProper& newStep(const StepProper& step);
    StepProper& newStep(StepProper&& step={});
    std::vector<StepProper>& newSteps(const std::vector<StepProper> &v);
    StepProper& getStep(size_t idx);
    const StepProper& getStep(size_t idx) const;
    std::vector<StepProper>& getSteps(void) noexcept;
    const std::vector<StepProper>& getSteps(void) const noexcept;
    size_t getNstep(void) const noexcept;

    void setName(const std::string &s);
    const std::string& getName(void) const noexcept;

    void setKPoints(const KPoints &k);
    KPoints& getKPoints(void) noexcept;
    const KPoints& getKPoints(void) const noexcept;
private:
    std::vector<StepProper> steps;
    std::string name;
    KPoints kpoints;
};
}
#endif // MOLECULE_H
