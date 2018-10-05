#ifndef MOLECULE_H
#define MOLECULE_H

#include "step.h"
#include "kpoints.h"
#include <list>

namespace Vipster {

class Molecule
{
public:
    Molecule(const std::string &name="New Molecule",unsigned long s=1);
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
