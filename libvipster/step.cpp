#include <global.h>
#include <step.h>
//TODO: remove debug headers
#include <iostream>

using namespace Vipster;

Step::Step():
    pse{std::make_shared<PseMap>()}
{
}

Step::Step(const std::shared_ptr<PseMap> &pse):
    pse{pse}
{
}

std::ostream& Vipster::operator<< (std::ostream& s, const Step& st)
{
    const Mat &v = st.cellvec;
    s << "Step:\n Atoms: " << st.getNat()
      <<"\n Types: " << st.getNtyp()
      << "\n Cell dimension: " << st.getCellDim()
      << "\n Vectors:\n [[" << v[0][0] << ", " << v[0][1] << ", " << v[0][2]
      << "]\n  [" << v[1][0] << ", " << v[1][1] << ", " << v[1][2]
      << "]\n  [" << v[2][0] << ", " << v[2][1] << ", " << v[2][2] << "]]"
      << "\n Comment: " << st.comment;
    return s;
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
    atoms.push_back(formatAtom(at,fmt,format));
    bonds_outdated = true;
}

void Step::newAtoms(const std::vector<Atom> &v)
{
    atoms.reserve(atoms.size()+v.size());
    atoms.insert(atoms.end(),v.begin(),v.end());
    bonds_outdated = true;
}

void Step::delAtom(size_t idx)
{
    if(idx>atoms.size()-1)throw std::out_of_range("Step::delAtom() : index is out of range");
    atoms.erase(atoms.begin()+idx);
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
    atoms.at(idx) = formatAtom(at,fmt,format);
    bonds_outdated = true;
}

const Atom& Step::getAtom(size_t idx) const
{
    return atoms.at(idx);
}

Atom Step::getAtom(size_t idx, AtomFmt fmt) const
{
    return formatAtom(atoms.at(idx),format,fmt);
}

const std::vector<Atom>& Step::getAtoms() const
{
    return atoms;
}

std::vector<Atom> Step::getAtoms(AtomFmt fmt) const
{
    return formatAtoms(atoms,format,fmt);
}

size_t Step::getNat() const noexcept
{
    return atoms.size();
}

AtomFmt Step::getFmt() const noexcept
{
    return format;
}

void Step::setFmt(AtomFmt fmt, bool scale)
{
    if(scale)
        atoms = formatAtoms(atoms,format,fmt);
    format = fmt;
}

Atom Step::formatAtom(Atom at, AtomFmt source, AtomFmt target) const
{
    if (source == target) return at;
    switch(source){
    case AtomFmt::Bohr:
        break;
    case AtomFmt::Angstrom:
        at.coord*=invbohr;
        break;
    case AtomFmt::Crystal:
        at.coord = at.coord * cellvec;
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
        at.coord = at.coord * invvec;
    case AtomFmt::Alat:
        at.coord /= celldim;
        break;
    }
    return at;
}

std::vector<Atom> Step::formatAtoms(std::vector<Atom> atvec, AtomFmt source, AtomFmt target) const
{
    if(source==target) return atvec;
    switch(source)
    {
    case AtomFmt::Bohr:
        break;
    case AtomFmt::Angstrom:
        for(Atom& at:atvec) { at.coord *= invbohr; }
        break;
    case AtomFmt::Crystal:
        for(Atom& at:atvec) { at.coord = at.coord * cellvec; }
    case AtomFmt::Alat:
        for(Atom& at:atvec) { at.coord *= celldim; }
        break;
    }
    switch(target)
    {
    case AtomFmt::Bohr:
        break;
    case AtomFmt::Angstrom:
        for(Atom& at:atvec) { at.coord *= bohrrad; }
        break;
    case AtomFmt::Crystal:
        for(Atom& at:atvec) { at.coord = at.coord * invvec; }
    case AtomFmt::Alat:
        for(Atom& at:atvec) { at.coord /= celldim; }
        break;
    }
    return atvec;
}

void Step::setCellDim(float cdm, bool scale)
{
    if(!(cdm>0))throw std::invalid_argument("Step::setCellDim() : cell-dimension needs to be positive");
    if(scale)
    {
        float ratio = cdm/celldim;
        for(Atom& at:atoms) { at.coord*=ratio; }
    }
    celldim = cdm;
    bonds_outdated = true;
}

void Step::setCellDim(float cdm, bool scale, AtomFmt fmt)
{
    if(fmt==AtomFmt::Angstrom){
        cdm *= invbohr;
    }
    setCellDim(cdm, scale);
}


float Step::getCellDim() const noexcept
{
    return celldim;
}

float Step::getCellDim(AtomFmt fmt) const noexcept
{
    if(fmt!=format){
        if(fmt==AtomFmt::Angstrom){
            return celldim * bohrrad;
        }else{
            return celldim * invbohr;
        }
    }else{
        return celldim;
    }
}

void Step::setCellVec(const Mat &mat, bool scale)
{
    Mat inv = Mat_inv(mat);
    std::vector<Atom> tatoms;
    if(scale){
        tatoms=formatAtoms(atoms,format,AtomFmt::Crystal);
    }
    cellvec = mat;
    invvec.swap(inv);
    if(scale){
        atoms=formatAtoms(tatoms,AtomFmt::Crystal,format);
    }
    bonds_outdated = true;
}

const Mat& Step::getCellVec() const noexcept
{
    return cellvec;
}

Vec Step::getCenter(bool com) const noexcept
{
    if(com && getNat()){
        Vec min{{std::numeric_limits<float>::max(),
                 std::numeric_limits<float>::max(),
                 std::numeric_limits<float>::max()}};
        Vec max{{std::numeric_limits<float>::min(),
                 std::numeric_limits<float>::min(),
                 std::numeric_limits<float>::min()}};
        for(const Atom& at:atoms){
            min[0]=std::min(min[0],at.coord[0]);
            min[1]=std::min(min[1],at.coord[1]);
            min[2]=std::min(min[2],at.coord[2]);
            max[0]=std::max(max[0],at.coord[0]);
            max[1]=std::max(max[1],at.coord[1]);
            max[2]=std::max(max[2],at.coord[2]);
        }
        return (min+max)/2;
    }else{
        return (cellvec[0]+cellvec[1]+cellvec[2])*celldim/2;
    }
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

void  Step::setComment(const std::string &s)
{
    comment = s;
}

const std::string& Step::getComment() const noexcept
{
    return comment;
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

size_t Step::getNbond() const noexcept
{
    return bonds.size();
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

void Step::checkBond(std::size_t i, std::size_t j, float cutfac, Vec dist, std::array<int, 3> offset) const
{
    if (dist[0] > cutfac) return;
    if (dist[1] > cutfac) return;
    if (dist[2] > cutfac) return;
    float dist_n = Vec_dot(dist, dist);
    if((0.57 < dist_n) && (dist_n < cutfac*cutfac)){
        bonds.push_back({i, j, std::sqrt(dist_n), offset[0], offset[1], offset[2]});
    }
}

void Step::setBondsCell(float cutfac) const
{
    bonds.clear();
    Vec x = cellvec[0] * celldim;
    Vec y = cellvec[1] * celldim;
    Vec z = cellvec[2] * celldim;
    Vec xy   = x+y;
    Vec xmy  = x-y;
    Vec xz   = x+z;
    Vec xmz  = x-z;
    Vec yz   = y+z;
    Vec ymz  = y-z;
    Vec xyz  = xy + z;
    Vec xymz = xy - z;
    Vec xmyz = xz - y;
    Vec mxyz = yz - x;
    std::array<int, 3> diff_v, crit_v;
    for(std::vector<Atom>::size_type i = 0; i != atoms.size(); ++i){
        float cut_i = (*pse)[atoms[i].name].bondcut;
        if (!cut_i) continue;
        for(std::vector<Atom>::size_type j = i+1; j != atoms.size(); ++j){
            float cut_j = (*pse)[atoms[j].name].bondcut;
            if (!cut_j) continue;
            float effcut = (cut_i + cut_j) * cutfac;
            Vec dist_v = atoms[i].coord - atoms[j].coord;
            if(format != AtomFmt::Crystal){
                dist_v = dist_v * invvec / celldim;
            }
            // TODO TODO TODO: vorzeichenfehler?! Bindungen werden teils falsch angezeigt!
            diff_v[0] = std::copysign(std::floor(std::abs(dist_v[0])), dist_v[0]);
            diff_v[1] = std::copysign(std::floor(std::abs(dist_v[1])), dist_v[1]);
            diff_v[2] = std::copysign(std::floor(std::abs(dist_v[2])), dist_v[2]);
            dist_v[0] = std::fmod(dist_v[0], 1);
            dist_v[1] = std::fmod(dist_v[1], 1);
            dist_v[2] = std::fmod(dist_v[2], 1);
            if(std::abs(dist_v[0]) < std::numeric_limits<float>::epsilon()){
                crit_v[0] = 0;
            }else{
                crit_v[0] = dist_v[0] < 0 ? -1 : 1;
            }
            if(std::abs(dist_v[1]) < std::numeric_limits<float>::epsilon()){
                crit_v[1] = 0;
            }else{
                crit_v[1] = dist_v[1] < 0 ? -1 : 1;
            }
            if(std::abs(dist_v[2]) < std::numeric_limits<float>::epsilon()){
                crit_v[2] = 0;
            }else{
                crit_v[2] = dist_v[2] < 0 ? -1 : 1;
            }
            if(!(crit_v[0]||crit_v[1]||crit_v[2])){
                // TODO: fail here? set flag? overlapping atoms!
                continue;
            }
            dist_v = dist_v * cellvec * celldim;
            // 0-vector
            checkBond(i, j, effcut, dist_v, diff_v);
            if(crit_v[0]){
                // x, -x
                checkBond(i, j, effcut, dist_v-crit_v[0]*x,
                          {{diff_v[0]+crit_v[0],diff_v[1],diff_v[2]}});
            }
            if(crit_v[1]){
                // y, -y
                checkBond(i, j, effcut, dist_v-crit_v[1]*y,
                          {{diff_v[0],diff_v[1]+crit_v[1],diff_v[2]}});
                if(crit_v[0]){
                    if(crit_v[0] == crit_v[1]){
                        // x+y, -x-y
                        checkBond(i, j, effcut, dist_v-crit_v[0]*xy,
                                  {{diff_v[0]+crit_v[0],
                                    diff_v[1]+crit_v[1],
                                    diff_v[2]}});
                    }else{
                        // x-y, -x+y
                        checkBond(i, j, effcut, dist_v-crit_v[0]*xmy,
                                  {{diff_v[0]+crit_v[0],
                                    diff_v[1]+crit_v[1],
                                    diff_v[2]}});
                    }
                }
            }
            if(crit_v[2]){
                // z, -z
                checkBond(i, j, effcut, dist_v-crit_v[2]*z,
                          {{diff_v[0],diff_v[1],diff_v[2]+crit_v[2]}});
                if(crit_v[0]){
                    if(crit_v[0] == crit_v[2]){
                        // x+z, -x-z
                        checkBond(i, j, effcut, dist_v-crit_v[0]*xz,
                                  {{diff_v[0]+crit_v[0],
                                    diff_v[1],
                                    diff_v[2]+crit_v[2]}});
                    }else{
                        // x-z, -x+z
                        checkBond(i, j, effcut, dist_v-crit_v[0]*xmz,
                                  {{diff_v[0]+crit_v[0],
                                    diff_v[1],
                                    diff_v[2]+crit_v[2]}});
                    }
                }
                if(crit_v[1]){
                    if(crit_v[1] == crit_v[2]){
                        // y+z, -y-z
                        checkBond(i, j, effcut, dist_v-crit_v[1]*yz,
                                  {{diff_v[0],
                                    diff_v[1]+crit_v[1],
                                    diff_v[2]+crit_v[2]}});
                    }else{
                        // y-z, -y+z
                        checkBond(i, j, effcut, dist_v-crit_v[1]*ymz,
                                  {{diff_v[0],
                                    diff_v[1]+crit_v[1],
                                    diff_v[2]+crit_v[2]}});
                    }
                    if(crit_v[0]){
                        if(crit_v[0] == crit_v[1]){
                            if(crit_v[0] == crit_v[2]){
                                // x+y+z, -x-y-z
                                checkBond(i, j, effcut, dist_v-crit_v[0]*xyz,
                                          {{diff_v[0]+crit_v[0],
                                            diff_v[1]+crit_v[1],
                                            diff_v[2]+crit_v[2]}});
                            }else{
                                // x+y-z, -x-y+z
                                checkBond(i, j, effcut, dist_v-crit_v[0]*xymz,
                                          {{diff_v[0]+crit_v[0],
                                            diff_v[1]+crit_v[1],
                                            diff_v[2]+crit_v[2]}});
                            }
                        }else{
                            if(crit_v[0] == crit_v[2]){
                                // x-y+z, -x+y-z
                                checkBond(i, j, effcut, dist_v-crit_v[0]*xmyz,
                                          {{diff_v[0]+crit_v[0],
                                            diff_v[1]+crit_v[1],
                                            diff_v[2]+crit_v[2]}});
                            }else{
                                // x-y-z, -x+y+z
                                checkBond(i, j, effcut, dist_v-crit_v[1]*mxyz,
                                          {{diff_v[0]+crit_v[0],
                                            diff_v[1]+crit_v[1],
                                            diff_v[2]+crit_v[2]}});
                            }
                        }
                    }
                }
            }
//            Vec ctemp;
//            ctemp[0] = dist_v[0]*invvec[0][0] + dist_v[1]*invvec[1][0] + dist_v[2]*invvec[2][0];
//            ctemp[1] = dist_v[0]*invvec[0][1] + dist_v[1]*invvec[1][1] + dist_v[2]*invvec[2][1];
//            ctemp[2] = dist_v[0]*invvec[0][2] + dist_v[1]*invvec[1][2] + dist_v[2]*invvec[2][2];
//            dist_v.swap(ctemp);
//            dist_v /= celldim;
//            dist_v += 0.5;
//            int xdiff = std::floor(dist_v[0]);
//            int ydiff = std::floor(dist_v[1]);
//            int zdiff = std::floor(dist_v[2]);
//            dist_v[0] = std::abs(std::fmod(dist_v[0], 1));
//            dist_v[1] = std::abs(std::fmod(dist_v[1], 1));
//            dist_v[2] = std::abs(std::fmod(dist_v[2], 1));
//            dist_v -= 0.5;
//            ctemp[0] = dist_v[0]*cellvec[0][0] + dist_v[1]*cellvec[1][0] + dist_v[2]*cellvec[2][0];
//            ctemp[1] = dist_v[0]*cellvec[0][1] + dist_v[1]*cellvec[1][1] + dist_v[2]*cellvec[2][1];
//            ctemp[2] = dist_v[0]*cellvec[0][2] + dist_v[1]*cellvec[1][2] + dist_v[2]*cellvec[2][2];
//            dist_v.swap(ctemp);
//            dist_v *= celldim;
//            float dist_n = Vec_length(dist_v);
//            if ((0.57 < dist_n) && (dist_n < effcut)){
//                bonds.push_back({i, j, std::sqrt(dist_n), xdiff, ydiff, zdiff});
//            }
        }
    }
    bonds_outdated = false;
    bondcut_factor = cutfac;
    bonds_level = BondLevel::Cell;
}
