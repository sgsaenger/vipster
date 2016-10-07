#include "step.h"
#include <limits>
#include <iostream>

using namespace Vipster;

Step::Step(PseMap pse):
    pse(pse),
    celldim{1.},
    cellvec{{ {{1.,0.,0.}},{{0.,1.,0.}},{{0.,0.,1.}} }},
    invvec{{ {{1.,0.,0.}},{{0.,1.,0.}},{{0.,0.,1.}} }},
    bonds_outdated{true},
    bondcut_factor{1.1}
{
}

void Step::newAtom(std::string name, Vec coord, float charge, std::array<bool,3> fix, bool hidden, Fmt fmt)
{
    Atom at{name,coord,charge,fix,hidden};
    atoms.push_back(formatAtom(at,fmt,Fmt::Bohr));
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

void Step::newAtom(Atom at, Fmt fmt)
{
    atoms.push_back(formatAtom(at,fmt,Fmt::Bohr));
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
    if(idx>atoms.size()-1)throw std::out_of_range("Atomview::delAtom() : index is out of range");
    atoms.erase(atoms.begin()+idx);
    bonds_outdated = true;
}

void Step::setAtom(size_t idx, std::string name, Vec coord, float charge, std::array<bool,3> fix, bool hidden, Fmt fmt)
{
    Atom at{name,coord,charge,fix,hidden};
    atoms.at(idx) = formatAtom(at,fmt,Fmt::Bohr);
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

void Step::setAtom(size_t idx, Atom at, Fmt fmt)
{
    atoms.at(idx) = formatAtom(at,fmt,Fmt::Bohr);
    bonds_outdated = true;
}

const Atom& Step::getAtom(size_t idx) const
{
    return atoms.at(idx);
}

Atom Step::getAtomFmt(size_t idx, Fmt fmt)
{
    return formatAtom(atoms.at(idx),Fmt::Bohr,fmt);
}

const std::vector<Atom>& Step::getAtoms() const noexcept
{
    return atoms;
}

std::vector<Atom> Step::getAtomsFmt(Fmt fmt)
{
    return formatAtoms(atoms,Fmt::Bohr,fmt);
}

size_t Step::getNat() const noexcept
{
    return atoms.size();
}

Atom Step::formatAtom(Atom at, Fmt source, Fmt target)
{
    switch(source){
    case Fmt::Bohr:
        break;
    case Fmt::Angstrom:
        at.coord*=invbohr;
        break;
    case Fmt::Crystal:
        Vec ctemp;
        ctemp[0] = at.coord[0]*this->cellvec[0][0]+
                   at.coord[1]*this->cellvec[1][0]+
                   at.coord[2]*this->cellvec[2][0];
        ctemp[1] = at.coord[0]*this->cellvec[0][1]+
                   at.coord[1]*this->cellvec[1][1]+
                   at.coord[2]*this->cellvec[2][1];
        ctemp[2] = at.coord[0]*this->cellvec[0][2]+
                   at.coord[1]*this->cellvec[1][2]+
                   at.coord[2]*this->cellvec[2][2];
        at.coord.swap(ctemp);
    case Fmt::Alat:
        at.coord*=this->celldim;
        break;
    }
    switch(target){
    case Fmt::Bohr:
        break;
    case Fmt::Angstrom:
        at.coord*=bohrrad;
        break;
    case Fmt::Crystal:
        Vec ctemp;
        ctemp[0] = at.coord[0]*this->invvec[0][0]+
                   at.coord[1]*this->invvec[1][0]+
                   at.coord[2]*this->invvec[2][0];
        ctemp[1] = at.coord[0]*this->invvec[0][1]+
                   at.coord[1]*this->invvec[1][1]+
                   at.coord[2]*this->invvec[2][1];
        ctemp[2] = at.coord[0]*this->invvec[0][2]+
                   at.coord[1]*this->invvec[1][2]+
                   at.coord[2]*this->invvec[2][2];
        at.coord.swap(ctemp);
    case Fmt::Alat:
        at.coord/=this->celldim;
        break;
    }
    return at;
}

std::vector<Atom> Step::formatAtoms(std::vector<Atom> atoms, Fmt source, Fmt target)
{
    switch(source)
    {
    case Fmt::Bohr:
        break;
    case Fmt::Angstrom:
        for(Atom at:atoms) { at.coord*=invbohr; }
        break;
    case Fmt::Crystal:
        Vec ctemp;
        for(Atom at:atoms)
        {
            ctemp[0] = at.coord[0]*this->cellvec[0][0]+
                       at.coord[1]*this->cellvec[1][0]+
                       at.coord[2]*this->cellvec[2][0];
            ctemp[1] = at.coord[0]*this->cellvec[0][1]+
                       at.coord[1]*this->cellvec[1][1]+
                       at.coord[2]*this->cellvec[2][1];
            ctemp[2] = at.coord[0]*this->cellvec[0][2]+
                       at.coord[1]*this->cellvec[1][2]+
                       at.coord[2]*this->cellvec[2][2];
            at.coord.swap(ctemp);
        }
    case Fmt::Alat:
        for(Atom at:atoms) { at.coord*=this->celldim; }
    }
    switch(target)
    {
    case Fmt::Bohr:
        break;
    case Fmt::Angstrom:
        for(Atom at:atoms) { at.coord*=bohrrad; }
        break;
    case Fmt::Crystal:
        Vec ctemp;
        for(Atom at:atoms)
        {
            ctemp[0] = at.coord[0]*this->invvec[0][0]+
                       at.coord[1]*this->invvec[1][0]+
                       at.coord[2]*this->invvec[2][0];
            ctemp[1] = at.coord[0]*this->invvec[0][1]+
                       at.coord[1]*this->invvec[1][1]+
                       at.coord[2]*this->invvec[2][1];
            ctemp[2] = at.coord[0]*this->invvec[0][2]+
                       at.coord[1]*this->invvec[1][2]+
                       at.coord[2]*this->invvec[2][2];
            at.coord.swap(ctemp);
        }
    case Fmt::Alat:
        for(Atom at:atoms) { at.coord/=this->celldim; }
    }
    return atoms;
}

void Step::setCellDim(float cdm, bool scale)
{
    if(!(cdm>0))throw std::invalid_argument("Atomview::setCellDim() : cell-dimension needs to be positive");
    if(scale)
    {
        float ratio = cdm/celldim;
        for(Atom at:atoms) { at.coord*=ratio; }
    }
    celldim = cdm;
    bonds_outdated = true;
}

float Step::getCellDim() const noexcept
{
    return celldim;
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
                +vec[1][0]*(vec[1][2]*vec[2][0]-vec[1][0]*vec[2][2])
                +vec[2][0]*(vec[1][0]*vec[2][1]-vec[1][1]*vec[2][0]);
    if(std::abs(det) < std::numeric_limits<float>::epsilon())
    {
        throw std::invalid_argument("Atomview::setCellVec() : invalid cell-vectors (singular matrix)");
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
        tatoms=formatAtoms(atoms,Fmt::Bohr,Fmt::Crystal);
    }
    cellvec.swap(vec);
    invvec.swap(inv);
    if(scale){
        atoms=formatAtoms(tatoms,Fmt::Crystal,Fmt::Bohr);
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
        for(Atom at:atoms){
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
    for(Atom at:atoms) { types.insert(at.name); }
    return types;
}

size_t Step::getNtyp() const noexcept
{
    return getTypes().size();
}

const std::array<std::vector<Bond>,8>& Step::getBonds() const
{
    return getBonds(bondcut_factor);
}

const std::array<std::vector<Bond>,8>& Step::getBonds(float cutfac) const
{
    if(bonds_outdated or (cutfac!=bondcut_factor))
    {
        Step::setBonds(cutfac);
    }
    return bonds;
}

void Step::setBonds(float cutfac) const
{
    float  effcut,dist_n;
    Vec   dist_v;
    auto    offsets = getBondOffsets();
    for(uint dir=0;dir<8;++dir){
        bonds[dir].clear();
        for(auto off: offsets[dir]){
            for(std::vector<Atom>::size_type i = 0; i != atoms.size(); ++i){
                float cut_i = pse[atoms[i].name].bondcut;
                if(cut_i < std::numeric_limits<float>::epsilon()) continue;
                for(std::vector<Atom>::size_type j = 0; j != atoms.size(); ++j){
                    if( (j<i) && (dir == 0)) continue;
                    if( j == i) continue;
                    float cut_j = pse[atoms[j].name].bondcut;
                    if(cut_j < std::numeric_limits<float>::epsilon()) continue;
                    effcut = (cut_i+cut_j)*cutfac;
                    Vec pos_i = atoms[i].coord;
                    Vec pos_j = atoms[j].coord;
                    dist_v[0] = pos_i[0] + off[0][0] - pos_j[0] - off[1][0];
                    if(dist_v[0]>effcut)continue;
                    dist_v[1] = pos_i[1] + off[0][1] - pos_j[1] - off[1][1];
                    if(dist_v[1]>effcut)continue;
                    dist_v[2] = pos_i[2] + off[0][2] - pos_j[2] - off[1][2];
                    if(dist_v[2]>effcut)continue;
                    dist_n = dist_v[0]*dist_v[0]+dist_v[1]*dist_v[1]+dist_v[2]*dist_v[2];
                    if((0.57<dist_n)&&(dist_n<effcut*effcut)){
                        bonds[dir].push_back({i,j,std::sqrt(dist_n)});
                    }
                }
            }
        }
    }
    bonds_outdated = false;
    bondcut_factor = cutfac;
}

std::array<std::vector<std::array<Vec,2>>,8> Step::getBondOffsets() const
{
    Vec   n{0.};
    Vec   vec01{cellvec[0]+cellvec[1]};
    Vec   vec02{cellvec[0]+cellvec[2]};
    Vec   vec12{cellvec[1]+cellvec[2]};
    Vec   vec012{vec01+cellvec[2]};
    std::array<Vec,2> orig{n, n}, x{cellvec[0], n}, y{cellvec[1], n}, z{cellvec[2], n};
    std::array<Vec,2> xy{vec01, n}, xmy{cellvec[0], cellvec[1]};
    std::array<Vec,2> xz{vec02, n}, xmz{cellvec[0], cellvec[2]};
    std::array<Vec,2> yz{vec12, n}, ymz{cellvec[1], cellvec[2]};
    std::array<Vec,2> xyz{vec012,n},xymz{vec01,cellvec[2]},xmyz{vec02,cellvec[1]},mxyz{vec12,cellvec[0]};
    std::array<std::vector<std::array<Vec,2>>,8> off_list{{
        {{orig}},
        {{x}},{{y}}, {{xy,xmy}},
        {{z}}, {{xz,xmz}}, {{yz,ymz}},
        {{xyz,xymz,xmyz,mxyz}}
    }};
    return off_list;
}
