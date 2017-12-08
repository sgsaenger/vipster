#ifndef MOLECULE_H
#define MOLECULE_H

#include "step.h"
#include "config.h"
#include "kpoints.h"

namespace Vipster {
class Molecule
{
public:
    Molecule(std::string name="New Molecule",unsigned long s=1);
    std::shared_ptr<PseMap> pse = std::make_shared<PseMap>();

    void setCellDimAll(float cdm, bool scale=false, AtomFmt fmt=AtomFmt::Bohr);
    void setCellVecAll(const Mat &mat, bool scale=false);
    void setFmtAll(AtomFmt fmt, bool scale=false);

    Step& newStep(const Step& step);
    Step& newStep(Step&& step={});
    std::vector<Step>& newSteps(const std::vector<Step> &v);
    Step& getStep(size_t idx);
    const Step& getStep(size_t idx) const;
    std::vector<Step>& getSteps(void) noexcept;
    const std::vector<Step>& getSteps(void) const noexcept;
    size_t getNstep(void) const noexcept;

    void setName(const std::string &s);
    const std::string& getName(void) const noexcept;

    void setKPoints(const KPoints &k);
    KPoints& getKPoints(void) noexcept;
    const KPoints& getKPoints(void) const noexcept;
private:
    std::vector<Step> steps;
    std::string name;
    KPoints kpoints;
};
}
#endif // MOLECULE_H
