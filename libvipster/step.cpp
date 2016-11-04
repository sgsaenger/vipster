#include <step.h>
#include <limits>
#include <cmath>

using namespace Vipster;

Step::Step():
    pse{std::make_shared<PseMap>(&Vipster::pse)}
{
}

Step::Step(const std::shared_ptr<PseMap> &pse):
    pse{pse}
{
}

void Step::newAtom(std::string name, Vec coord, float charge, std::array<bool,3> fix, bool hidden, AtomFmt fmt)
{
    Atom at{name,coord,charge,fix,hidden};
    atoms.push_back(formatAtom(at,fmt,AtomFmt::Bohr));
    bonds_outdated = true;
}

void Step::newAtom(const Atom& at)
{
    atoms.push_back(at);
    bonds_outdated = true;
}

void Step::newAtom(Atom &&at)
{
    atoms.push_back(std::move(at));
    bonds_outdated = true;
}

void Step::newAtom(Atom at, AtomFmt fmt)
{
    atoms.push_back(formatAtom(at,fmt,AtomFmt::Bohr));
    bonds_outdated = true;
}

void Step::newAtoms(size_t count)
{
    Atom temp {"C",{0.,0.,0.},0.,{false,false,false},false};
    atoms.reserve(atoms.size()+count);
    for(size_t i=0;i!=count;i++)
    {
        atoms.push_back(temp);
    }
    bonds_outdated = true;
}

void Step::delAtom(size_t idx)
{
    if(idx>atoms.size()-1)throw std::out_of_range("Step::delAtom() : index is out of range");
    atoms.erase(atoms.begin()+idx);
    bonds_outdated = true;
}

void Step::setAtom(size_t idx, std::string name, Vec coord, float charge, std::array<bool,3> fix, bool hidden, AtomFmt fmt)
{
    Atom at{name,coord,charge,fix,hidden};
    atoms.at(idx) = formatAtom(at,fmt,AtomFmt::Bohr);
    bonds_outdated = true;
}

void Step::setAtom(size_t idx, const Atom &at)
{
    atoms.at(idx) = at;
    bonds_outdated = true;
}

void Step::setAtom(size_t idx, Atom &&at)
{
    atoms.at(idx) = std::move(at);
    bonds_outdated = true;
}

void Step::setAtom(size_t idx, Atom at, AtomFmt fmt)
{
    atoms.at(idx) = formatAtom(at,fmt,AtomFmt::Bohr);
    bonds_outdated = true;
}

const Atom& Step::getAtom(size_t idx) const
{
    return atoms.at(idx);
}

Atom Step::getAtomFmt(size_t idx, AtomFmt fmt)
{
    return formatAtom(atoms.at(idx),AtomFmt::Bohr,fmt);
}

const std::vector<Atom>& Step::getAtoms() const noexcept
{
    return atoms;
}

std::vector<Atom> Step::getAtomsFmt(AtomFmt fmt)
{
    return formatAtoms(atoms,AtomFmt::Bohr,fmt);
}

size_t Step::getNat() const noexcept
{
    return atoms.size();
}

Atom Step::formatAtom(Atom at, AtomFmt source, AtomFmt target)
{
    if (source == target) return at;
    switch(source){
    case AtomFmt::Bohr:
        break;
    case AtomFmt::Angstrom:
        at.coord*=invbohr;
        break;
    case AtomFmt::Crystal:
        Vec ctemp;
        ctemp[0] = at.coord[0]*cellvec[0][0]+
                   at.coord[1]*cellvec[1][0]+
                   at.coord[2]*cellvec[2][0];
        ctemp[1] = at.coord[0]*cellvec[0][1]+
                   at.coord[1]*cellvec[1][1]+
                   at.coord[2]*cellvec[2][1];
        ctemp[2] = at.coord[0]*cellvec[0][2]+
                   at.coord[1]*cellvec[1][2]+
                   at.coord[2]*cellvec[2][2];
        at.coord.swap(ctemp);
    case AtomFmt::Alat:
        at.coord *= celldim;
        break;
    }
    switch(target){
    case AtomFmt::Bohr:
        break;
    case AtomFmt::Angstrom:
        at.coord*=bohrrad;
        break;
    case AtomFmt::Crystal:
        Vec ctemp;
        ctemp[0] = at.coord[0]*invvec[0][0]+
                   at.coord[1]*invvec[1][0]+
                   at.coord[2]*invvec[2][0];
        ctemp[1] = at.coord[0]*invvec[0][1]+
                   at.coord[1]*invvec[1][1]+
                   at.coord[2]*invvec[2][1];
        ctemp[2] = at.coord[0]*invvec[0][2]+
                   at.coord[1]*invvec[1][2]+
                   at.coord[2]*invvec[2][2];
        at.coord.swap(ctemp);
    case AtomFmt::Alat:
        at.coord /= celldim;
        break;
    }
    return at;
}

std::vector<Atom> Step::formatAtoms(std::vector<Atom> atoms, AtomFmt source, AtomFmt target)
{
    switch(source)
    {
    case AtomFmt::Bohr:
        break;
    case AtomFmt::Angstrom:
        for(Atom& at:atoms) { at.coord*=invbohr; }
        break;
    case AtomFmt::Crystal:
        Vec ctemp;
        for(Atom& at:atoms)
        {
            ctemp[0] = at.coord[0]*cellvec[0][0]+
                       at.coord[1]*cellvec[1][0]+
                       at.coord[2]*cellvec[2][0];
            ctemp[1] = at.coord[0]*cellvec[0][1]+
                       at.coord[1]*cellvec[1][1]+
                       at.coord[2]*cellvec[2][1];
            ctemp[2] = at.coord[0]*cellvec[0][2]+
                       at.coord[1]*cellvec[1][2]+
                       at.coord[2]*cellvec[2][2];
            at.coord.swap(ctemp);
        }
    case AtomFmt::Alat:
        for(Atom& at:atoms) { at.coord*=celldim; }
        break;
    }
    switch(target)
    {
    case AtomFmt::Bohr:
        break;
    case AtomFmt::Angstrom:
        for(Atom& at:atoms) { at.coord*=bohrrad; }
        break;
    case AtomFmt::Crystal:
        Vec ctemp;
        for(Atom& at:atoms)
        {
            ctemp[0] = at.coord[0]*invvec[0][0]+
                       at.coord[1]*invvec[1][0]+
                       at.coord[2]*invvec[2][0];
            ctemp[1] = at.coord[0]*invvec[0][1]+
                       at.coord[1]*invvec[1][1]+
                       at.coord[2]*invvec[2][1];
            ctemp[2] = at.coord[0]*invvec[0][2]+
                       at.coord[1]*invvec[1][2]+
                       at.coord[2]*invvec[2][2];
            at.coord.swap(ctemp);
        }
    case AtomFmt::Alat:
        for(Atom& at:atoms) { at.coord/=celldim; }
        break;
    }
    return atoms;
}

void Step::setCellDim(float cdm, bool scale, AtomFmt fmt)
{
    if(fmt==AtomFmt::Angstrom){
        cdm *= invbohr;
    }
    if(!(cdm>0))throw std::invalid_argument("Step::setCellDim() : cell-dimension needs to be positive");
    if(scale)
    {
        float ratio = cdm/celldim;
        for(Atom& at:atoms) { at.coord*=ratio; }
    }
    celldim = cdm;
    bonds_outdated = true;
}

float Step::getCellDim(AtomFmt fmt) const noexcept
{
    if(fmt==AtomFmt::Angstrom){
        return celldim * bohrrad;
    }else{
        return celldim;
    }
}

void Step::setCellVec(float v11, float v12, float v13, float v21, float v22, float v23, float v31, float v32, float v33, bool scale)
{
    setCellVec(std::array<Vec,3>{{{{v11,v12,v13}},{{v21,v22,v23}},{{v31,v32,v33}}}},scale);
}

void Step::setCellVec(Vec v1, Vec v2, Vec v3, bool scale)
{
    setCellVec(std::array<Vec,3>{v1,v2,v3},scale);
}

void Step::setCellVec(std::array<Vec, 3> vec, bool scale)
{
    float det = vec[0][0]*(vec[1][1]*vec[2][2]-vec[1][2]*vec[2][1])
               +vec[0][1]*(vec[1][2]*vec[2][0]-vec[1][0]*vec[2][2])
               +vec[0][2]*(vec[1][0]*vec[2][1]-vec[1][1]*vec[2][0]);
    if(std::abs(det) < std::numeric_limits<float>::epsilon())
    {
        throw std::invalid_argument("Step::setCellVec() : invalid cell-vectors (singular matrix)");
    }
    float invdet = 1/det;
    std::array<Vec,3> inv;
    inv[0][0] = (vec[1][1]*vec[2][2]-vec[2][1]*vec[1][2])*invdet;
    inv[0][1] = (vec[0][2]*vec[2][1]-vec[0][1]*vec[2][2])*invdet;
    inv[0][2] = (vec[0][1]*vec[1][2]-vec[0][2]*vec[1][1])*invdet;
    inv[1][0] = (vec[1][2]*vec[2][0]-vec[1][0]*vec[2][2])*invdet;
    inv[1][1] = (vec[0][0]*vec[2][2]-vec[0][2]*vec[2][0])*invdet;
    inv[1][2] = (vec[1][0]*vec[0][2]-vec[0][0]*vec[1][2])*invdet;
    inv[2][0] = (vec[1][0]*vec[2][1]-vec[2][0]*vec[1][1])*invdet;
    inv[2][1] = (vec[2][0]*vec[0][1]-vec[0][0]*vec[2][1])*invdet;
    inv[2][2] = (vec[0][0]*vec[1][1]-vec[1][0]*vec[0][1])*invdet;
    std::vector<Atom> tatoms;
    if(scale){
        tatoms=formatAtoms(atoms,AtomFmt::Bohr,AtomFmt::Crystal);
    }
    cellvec.swap(vec);
    invvec.swap(inv);
    if(scale){
        atoms=formatAtoms(tatoms,AtomFmt::Crystal,AtomFmt::Bohr);
    }
    bonds_outdated = true;
}

const std::array<Vec,3>& Step::getCellVec() const noexcept
{
    return cellvec;
}

Vec Step::getCenter(bool com) const
{
    Vec temp{0.,0.,0.};
    if(com){
        Vec min{0.},max{0.};
        for(const Atom& at:atoms){
            min[0]=(min[0]<at.coord[0])?min[0]:at.coord[0];
            min[1]=(min[1]<at.coord[1])?min[1]:at.coord[1];
            min[2]=(min[2]<at.coord[2])?min[2]:at.coord[2];
            max[0]=(max[0]>at.coord[0])?max[0]:at.coord[0];
            max[1]=(max[1]>at.coord[1])?max[1]:at.coord[1];
            max[2]=(max[2]>at.coord[2])?max[2]:at.coord[2];
        }
        temp=(min+max)/2;
    }else{
        temp=(cellvec[0]+cellvec[1]+cellvec[2])*celldim/2;
    }
    return temp;
}

std::set<std::string> Step::getTypes() const noexcept
{
    std::set<std::string> types;
    for(const Atom& at:atoms) { types.insert(at.name); }
    return types;
}

size_t Step::getNtyp() const noexcept
{
    return getTypes().size();
}

const std::vector<Bond>& Step::getBonds() const
{
    return getBonds(bondcut_factor);
}

const std::vector<Bond>& Step::getBonds(float cutfac) const
{
    if(bonds_outdated or (cutfac!=bondcut_factor) or (bonds_level<BondLevel::Molecule))
    {
        Step::setBonds(cutfac);
    }
    return bonds;
}

const std::vector<Bond>& Step::getBondsCell() const
{
    return getBondsCell(bondcut_factor);
}

const std::vector<Bond>& Step::getBondsCell(float cutfac) const
{
    if(bonds_outdated or (cutfac!=bondcut_factor) or (bonds_level<BondLevel::Cell))
    {
        Step::setBondsCell(cutfac);
    }
    return bonds;
}

void Step::setBonds(float cutfac) const
{
    bonds.clear();
    for(std::vector<Atom>::size_type i = 0; i != atoms.size(); ++i){
        float cut_i = (*pse)[atoms[i].name].bondcut;
        if (!cut_i) continue;
        for(std::vector<Atom>::size_type j = i+1; j != atoms.size(); ++j){
            float cut_j = (*pse)[atoms[j].name].bondcut;
            if (!cut_j) continue;
            float effcut = (cut_i + cut_j) * cutfac;
            Vec pos_i = atoms[i].coord;
            Vec pos_j = atoms[j].coord;
            Vec dist_v;
            dist_v[0] = pos_i[0] - pos_j[0];
            if (dist_v[0] > effcut) continue;
            dist_v[1] = pos_i[1] - pos_j[1];
            if (dist_v[1] > effcut) continue;
            dist_v[2] = pos_i[2] - pos_j[2];
            if (dist_v[2] > effcut) continue;
            float dist_n = Vec_dot(dist_v, dist_v);
            if((0.57 < dist_n) && (dist_n < effcut*effcut)){
                bonds.push_back({i, j, std::sqrt(dist_n), 0, 0, 0});
            }
        }
    }
    bonds_outdated = false;
    bondcut_factor = cutfac;
    bonds_level = BondLevel::Molecule;
}

void Step::setBondsCell(float cutfac) const
{
    bonds.clear();
    for(std::vector<Atom>::size_type i = 0; i != atoms.size(); ++i){
        float cut_i = (*pse)[atoms[i].name].bondcut;
        if (!cut_i) continue;
        for(std::vector<Atom>::size_type j = i+1; j != atoms.size(); ++j){
            float cut_j = (*pse)[atoms[j].name].bondcut;
            if (!cut_j) continue;
            float effcut = (cut_i + cut_j) * cutfac;
            Vec dist_v = atoms[i].coord - atoms[j].coord;
            Vec ctemp;
            ctemp[0] = dist_v[0]*invvec[0][0] + dist_v[1]*invvec[1][0] + dist_v[2]*invvec[2][0];
            ctemp[1] = dist_v[0]*invvec[0][1] + dist_v[1]*invvec[1][1] + dist_v[2]*invvec[2][1];
            ctemp[2] = dist_v[0]*invvec[0][2] + dist_v[1]*invvec[1][2] + dist_v[2]*invvec[2][2];
            dist_v.swap(ctemp);
            dist_v /= celldim;
            dist_v += 0.5;
            int xdiff = std::floor(dist_v[0]);
            int ydiff = std::floor(dist_v[1]);
            int zdiff = std::floor(dist_v[2]);
            dist_v[0] = std::abs(std::fmod(dist_v[0], 1));
            dist_v[1] = std::abs(std::fmod(dist_v[1], 1));
            dist_v[2] = std::abs(std::fmod(dist_v[2], 1));
            dist_v -= 0.5;
            ctemp[0] = dist_v[0]*cellvec[0][0] + dist_v[1]*cellvec[1][0] + dist_v[2]*cellvec[2][0];
            ctemp[1] = dist_v[0]*cellvec[0][1] + dist_v[1]*cellvec[1][1] + dist_v[2]*cellvec[2][1];
            ctemp[2] = dist_v[0]*cellvec[0][2] + dist_v[1]*cellvec[1][2] + dist_v[2]*cellvec[2][2];
            dist_v.swap(ctemp);
            dist_v *= celldim;
            float dist_n = Vec_length(dist_v);
            if ((0.57 < dist_n) && (dist_n < effcut)){
                bonds.push_back({i, j, std::sqrt(dist_n), xdiff, ydiff, zdiff});
            }
        }
    }
    bonds_outdated = false;
    bondcut_factor = cutfac;
    bonds_level = BondLevel::Cell;
}
