#include "molecule.h"

using namespace Vipster;


Molecule::Molecule(std::string name, ulong s):
    name{name},
    kpoints{}
{
    for(std::vector<Step>::size_type i=0;i!=s;++i){
        steps.emplace_back(pse);
    }
}

std::ostream& Vipster::operator<< (std::ostream& s, const Molecule& m)
{
    s << "Molecule:\n Name: " << m.getName() << "\n Steps: " << m.getNstep()
      << "\n Active K-Point: " << m.getKPoints();
    return s;
}


void Molecule::setCellDimAll(float cdm, bool scale, AtomFmt fmt)
{
    for(Step& s:steps){
        s.setCellDim(cdm,scale,fmt);
    }
}

void Molecule::setCellVecAll(const Mat &mat, bool scale)
{
    Mat inv = Mat_inv(mat);
    if(scale){
        std::vector<Atom> tatoms;
        for(Step& s:steps){
            tatoms = s.formatAtoms(s.atoms, AtomFmt::Bohr, AtomFmt::Crystal);
            s.cellvec = mat;
            s.invvec = inv;
            s.atoms = s.formatAtoms(tatoms, AtomFmt::Crystal, AtomFmt::Bohr);
            s.bonds_outdated = true;
        }
    }else{
        for(Step& s:steps){
            s.cellvec = mat;
            s.invvec = inv;
        }
    }
}

void Molecule::setFmtAll(AtomFmt fmt, bool scale)
{
    for(Step& s:steps){
        s.setFmt(fmt, scale);
    }
}

void Molecule::newStep(const Step &step)
{
    steps.push_back(step);
    steps.back().pse = pse;
}

void Molecule::newStep(Step &&step)
{
    steps.push_back(std::move(step));
    steps.back().pse = pse;
}

void Molecule::newSteps(const std::vector<Step> &v)
{
    steps.reserve(steps.size()+v.size());
    std::vector<Step>::iterator pos = steps.end();
    steps.insert(steps.end(),v.begin(),v.end());
    for(;pos!=steps.end();++pos)
    {
        pos->pse = pse;
    }
}

Step& Molecule::getStep(size_t idx)
{
    return steps.at(idx);
}

const Step& Molecule::getStep(size_t idx) const
{
    return steps.at(idx);
}

std::vector<Step>& Molecule::getSteps() noexcept
{
    return steps;
}

const std::vector<Step>& Molecule::getSteps() const noexcept
{
    return steps;
}

size_t Molecule::getNstep() const noexcept
{
    return steps.size();
}

void Molecule::setName(const std::string &s)
{
    name = s;
}

const std::string& Molecule::getName(void)const noexcept
{
    return name;
}

void Molecule::setKPoints(const KPoints &k)
{
    kpoints = k;
}

const KPoints& Molecule::getKPoints() const noexcept
{
    return kpoints;
}

KPoints& Molecule::getKPoints() noexcept
{
    return kpoints;
}
