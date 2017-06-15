#ifndef MOLECULE_H
#define MOLECULE_H

#include <step.h>
#include <config.h>
#include <kpoints.h>

namespace Vipster {
class Molecule
{
public:
    Molecule(std::string name="New Molecule",ulong s=1);
    std::shared_ptr<PseMap> pse = std::make_shared<PseMap>();
    friend std::ostream& operator<< (std::ostream& s, const Molecule& st);

    void setCellDimAll(float cdm, bool scale=false, AtomFmt fmt=AtomFmt::Bohr);

    void newStep(const Step& step);
    void newStep(Step&& step={});
    void newSteps(const std::vector<Step> &v);
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
std::ostream& operator<< (std::ostream& s, const Molecule& st);
}
#endif // MOLECULE_H
