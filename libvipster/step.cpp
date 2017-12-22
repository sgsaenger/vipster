#include "step.h"

using namespace Vipster;

Step::Step(const std::shared_ptr<PseMap> &pse,
           const std::shared_ptr<std::string> &comment,
           const std::shared_ptr<std::vector<std::string> > &at_name,
           const std::shared_ptr<std::vector<Vec> > &at_coord,
           const std::shared_ptr<std::vector<float> > &at_charge,
           const std::shared_ptr<std::vector<FixVec> > &at_fix,
           const std::shared_ptr<std::vector<char> > &at_hidden,
           const std::shared_ptr<float> &celldim,
           const std::shared_ptr<Mat> &cellvec,
           const std::shared_ptr<Mat> &invvec,
           const std::shared_ptr<BondLevel> &bonds_level,
           const std::shared_ptr<float> &bondcut_factor,
           const std::shared_ptr<std::vector<Bond> > &bonds)
    : pse{pse}, comment{comment}, at_name{at_name}, at_coord{at_coord}, at_charge{at_charge},
      at_fix{at_fix}, at_hidden{at_hidden}, celldim{celldim}, cellvec{cellvec}, invvec{invvec},
      bonds_level{bonds_level}, bondcut_factor{bondcut_factor}, bonds{bonds} {}

StepProper::StepProper()
    : StepProper(std::make_shared<PseMap>()) {}

StepProper::StepProper(const std::shared_ptr<PseMap> &pse)
    : Step{pse, std::make_shared<std::string>(),
           std::make_shared<std::vector<std::string>>(),
           std::make_shared<std::vector<Vec>>(),
           std::make_shared<std::vector<float>>(),
           std::make_shared<std::vector<FixVec>>(),
           std::make_shared<std::vector<char>>(),
           std::make_shared<float>(1),
           std::make_shared<Mat>(Mat{{{{1,0,0}},{{0,1,0}},{{0,0,1}}}}),
           std::make_shared<Mat>(Mat{{{{1,0,0}},{{0,1,0}},{{0,0,1}}}}),
           std::make_shared<BondLevel>(BondLevel::None),
           std::make_shared<float>(1.1),
           std::make_shared<std::vector<Bond>>()} {}

//TODO:
//enable formatting!

//Step::Step(const Step &s)
//    : asAlat{this, AtomFmt::Alat}, asAngstrom{this, AtomFmt::Angstrom},
//      asBohr{this, AtomFmt::Bohr}, asCrystal{this, AtomFmt::Crystal},
//      pse{s.pse} {}

//Step& Step::operator=(const Step& s)
//{
//    pse = s.pse;
//    return *this;
//}

//Step::StepFormatter::StepFormatter(Step *step, AtomFmt fmt)
//    :step{step}, fmt{fmt} {}

//AtomProper Step::StepFormatter::operator [](size_t idx)
//{
//    //TODO: formatting =D
//    return {step->at_name[idx], step->at_coord[idx], step->at_charge[idx],
//            step->at_fix[idx], step->at_hidden[idx]};
//}

Step::iterator::iterator(const Step *step, size_t idx)
    :step{const_cast<Step*>(step)}, idx{idx}, at{(*step)[idx]} {}

Step::iterator& Step::iterator::operator++()
{
    at = (*step)[++idx];
    return *this;
}

Step::iterator Step::iterator::operator++(int)
{
    Step::iterator tmp(*this);
    operator++();
    return tmp;
}

Atom& Step::iterator::operator*()
{
    return at;
}

bool Step::iterator::operator!=(const Step::iterator& rhs)
{
    return ((step != rhs.step) || (idx != rhs.idx));
}

Step::iterator Step::begin()
{
    return Step::iterator(this, 0);
}

const Step::iterator Step::begin() const
{
    return Step::iterator(this, 0);
}

Step::iterator Step::end()
{
    return Step::iterator(this, getNat());
}

const Step::iterator Step::end() const
{
    return Step::iterator(this, getNat());
}

void Step::newAtom(const Atom &at)
{
    at_name->push_back(at.name);
    at_coord->push_back(at.coord);
    at_charge->push_back(at.charge);
    at_fix->push_back(at.fix);
    at_hidden->push_back(at.hidden);
    bonds_outdated = true;
}

void Step::newAtoms(size_t i)
{
    size_t oldNat = getNat();
    size_t nat = oldNat + i;
    at_name->resize(nat);
    at_coord->resize(nat);
    at_charge->resize(nat);
    at_fix->resize(nat);
    at_hidden->resize(nat);
    bonds_outdated = true;
}

void Step::delAtom(size_t idx)
{
    at_name->erase(at_name->begin()+idx);
    at_coord->erase(at_coord->begin()+idx);
    at_charge->erase(at_charge->begin()+idx);
    at_fix->erase(at_fix->begin()+idx);
    at_hidden->erase(at_hidden->begin()+idx);
    bonds_outdated = true;
}

Atom Step::operator[](size_t idx)
{
    return {&(*at_name)[idx], &(*at_coord)[idx], &(*at_charge)[idx],
            &(*at_fix)[idx], &(*at_hidden)[idx], &at_changed};
}

const Atom Step::operator[](size_t idx) const
{
    return {&(*at_name)[idx], &(*at_coord)[idx], &(*at_charge)[idx],
            &(*at_fix)[idx], &(*at_hidden)[idx], &at_changed};
}

//AtomProper Step::getAtomFmt(size_t idx, AtomFmt fmt) const
//{
//    return formatAtom(atoms.at(idx),format,fmt);
//}

size_t Step::getNat() const noexcept
{
    return at_coord->size();
}

AtomFmt Step::getFmt() const noexcept
{
    return fmt;
}

void StepProper::setFmt(AtomFmt format, bool scale)
{
    //TODO
//    if(scale)
//        atoms = formatAtoms(atoms,format,fmt);
    fmt = format;
}

//AtomProper Step::formatAtom(AtomProper at, AtomFmt source, AtomFmt target) const
//{
//    if (source == target) return at;
//    switch(source){
//    case AtomFmt::Bohr:
//        break;
//    case AtomFmt::Angstrom:
//        at.coord = at.coord * invbohr;
//        break;
//    case AtomFmt::Crystal:
//        at.coord = at.coord * cellvec;
//        [[fallthrough]];
//    case AtomFmt::Alat:
//        at.coord = at.coord * celldim;
//        break;
//    }
//    switch(target){
//    case AtomFmt::Bohr:
//        break;
//    case AtomFmt::Angstrom:
//        at.coord = at.coord * bohrrad;
//        break;
//    case AtomFmt::Crystal:
//        at.coord = at.coord * invvec;
//        [[fallthrough]];
//    case AtomFmt::Alat:
//        at.coord = at.coord / celldim;
//        break;
//    }
//    return at;
//}

//std::vector<AtomProper> Step::formatAtoms(std::vector<AtomProper> atvec, AtomFmt source, AtomFmt target) const
//{
//    if(source==target) return atvec;
//    switch(source)
//    {
//    case AtomFmt::Bohr:
//        break;
//    case AtomFmt::Angstrom:
//        for(AtomProper& at:atvec) { at.coord *= invbohr; }
//        break;
//    case AtomFmt::Crystal:
//        for(AtomProper& at:atvec) { at.coord = at.coord * cellvec; }
//        [[fallthrough]];
//    case AtomFmt::Alat:
//        for(AtomProper& at:atvec) { at.coord *= celldim; }
//        break;
//    }
//    switch(target)
//    {
//    case AtomFmt::Bohr:
//        break;
//    case AtomFmt::Angstrom:
//        for(AtomProper& at:atvec) { at.coord *= bohrrad; }
//        break;
//    case AtomFmt::Crystal:
//        for(AtomProper& at:atvec) { at.coord = at.coord * invvec; }
//        [[fallthrough]];
//    case AtomFmt::Alat:
//        for(AtomProper& at:atvec) { at.coord /= celldim; }
//        break;
//    }
//    return atvec;
//}

void StepProper::setCellDim(float cdm, bool scale)
{
    if(!(cdm>0))throw Error("Step::setCellDim() : cell-dimension needs to be positive");
    if(scale)
    {
        float ratio = cdm / *celldim;
        for(auto& at:*this) {at.coord *= ratio;}
    }
    *celldim = cdm;
    bonds_outdated = true;
}

float StepProper::getCellDim() const noexcept
{
    return *celldim;
}

void Step::setCellVec(const Mat &mat, bool scale)
{
    Mat inv = Mat_inv(mat);
    //TODO
//    if(scale && (format!=AtomFmt::Crystal)){
//        atoms = formatAtoms(atoms,format,AtomFmt::Crystal);
//    }
    *cellvec = mat;
    invvec->swap(inv);
//    if(scale && (format!=AtomFmt::Crystal)){
//        atoms = formatAtoms(atoms,AtomFmt::Crystal,format);
//    }
    bonds_outdated = true;
}

Mat Step::getCellVec() const noexcept
{
    return *cellvec;
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
        for(auto& at:*this){
            min[0]=std::min(min[0],at.coord[0]);
            min[1]=std::min(min[1],at.coord[1]);
            min[2]=std::min(min[2],at.coord[2]);
            max[0]=std::max(max[0],at.coord[0]);
            max[1]=std::max(max[1],at.coord[1]);
            max[2]=std::max(max[2],at.coord[2]);
        }
        return (min+max)/2;
    }else{
        return ((*cellvec)[0]+(*cellvec)[1]+(*cellvec)[2]) * *celldim / 2;
    }
    //TODO: format return value!
}

std::set<std::string> Step::getTypes() const noexcept
{
    std::set<std::string> types;
    for(const auto& at:*this) { types.insert(at.name); }
    return types;
}

size_t Step::getNtyp() const noexcept
{
    return getTypes().size();
}

const std::vector<Bond>& Step::getBonds() const
{
    return getBonds(*bondcut_factor);
}

const std::vector<Bond>& Step::getBonds(float cutfac) const
{
    //TODO
//    if(bonds_outdated or (cutfac!=bondcut_factor) or (bonds_level<BondLevel::Molecule))
//    {
//        Step::setBonds(cutfac);
//    }
    return *bonds;
}

const std::vector<Bond>& Step::getBondsCell() const
{
    return getBondsCell(*bondcut_factor);
}

const std::vector<Bond>& Step::getBondsCell(float cutfac) const
{
    //TODO
//    if(bonds_outdated or (cutfac!=bondcut_factor) or (bonds_level<BondLevel::Cell))
//    {
//        Step::setBondsCell(cutfac);
//    }
    return *bonds;
}

size_t Step::getNbond() const noexcept
{
    return bonds->size();
}

//void Step::setBonds(float cutfac) const
//{
//    bonds.clear();
//    for(std::vector<AtomProper>::size_type i = 0; i != atoms.size(); ++i){
//        float cut_i = (*pse)[atoms[i].name].bondcut;
//        if (!cut_i) continue;
//        for(std::vector<AtomProper>::size_type j = i+1; j != atoms.size(); ++j){
//            float cut_j = (*pse)[atoms[j].name].bondcut;
//            if (!cut_j) continue;
//            float effcut = (cut_i + cut_j) * cutfac;
//            Vec pos_i = atoms[i].coord;
//            Vec pos_j = atoms[j].coord;
//            Vec dist_v;
//            dist_v[0] = pos_i[0] - pos_j[0];
//            if (dist_v[0] > effcut) continue;
//            dist_v[1] = pos_i[1] - pos_j[1];
//            if (dist_v[1] > effcut) continue;
//            dist_v[2] = pos_i[2] - pos_j[2];
//            if (dist_v[2] > effcut) continue;
//            float dist_n = Vec_dot(dist_v, dist_v);
//            if((0.57 < dist_n) && (dist_n < effcut*effcut)){
//                bonds.push_back({i, j, std::sqrt(dist_n), 0, 0, 0});
//            }
//        }
//    }
//    bonds_outdated = false;
//    bondcut_factor = cutfac;
//    bonds_level = BondLevel::Molecule;
//}

void Step::checkBond(std::size_t i, std::size_t j, float cutfac, Vec dist, std::array<int, 3> offset) const
{
    if (dist[0] > cutfac) return;
    if (dist[1] > cutfac) return;
    if (dist[2] > cutfac) return;
    float dist_n = Vec_dot(dist, dist);
    if((0.57 < dist_n) && (dist_n < cutfac*cutfac)){
        bonds->push_back({i, j, std::sqrt(dist_n), offset[0], offset[1], offset[2]});
    }
}

//void Step::setBondsCell(float cutfac) const
//{
//    bonds.clear();
//    Vec x = cellvec[0] * celldim;
//    Vec y = cellvec[1] * celldim;
//    Vec z = cellvec[2] * celldim;
//    Vec xy   = x+y;
//    Vec xmy  = x-y;
//    Vec xz   = x+z;
//    Vec xmz  = x-z;
//    Vec yz   = y+z;
//    Vec ymz  = y-z;
//    Vec xyz  = xy + z;
//    Vec xymz = xy - z;
//    Vec xmyz = xz - y;
//    Vec mxyz = yz - x;
//    std::array<int, 3> diff_v, crit_v;
//    for(std::vector<AtomProper>::size_type i = 0; i != atoms.size(); ++i){
//        float cut_i = (*pse)[atoms[i].name].bondcut;
//        if (!cut_i) continue;
//        for(std::vector<AtomProper>::size_type j = i+1; j != atoms.size(); ++j){
//            float cut_j = (*pse)[atoms[j].name].bondcut;
//            if (!cut_j) continue;
//            float effcut = (cut_i + cut_j) * cutfac;
//            Vec dist_v = atoms[i].coord - atoms[j].coord;
//            if(format != AtomFmt::Crystal){
//                dist_v = dist_v * invvec / celldim;
//            }
//            diff_v[0] = std::copysign(std::floor(std::abs(dist_v[0])), dist_v[0]);
//            diff_v[1] = std::copysign(std::floor(std::abs(dist_v[1])), dist_v[1]);
//            diff_v[2] = std::copysign(std::floor(std::abs(dist_v[2])), dist_v[2]);
//            dist_v[0] = std::fmod(dist_v[0], 1);
//            dist_v[1] = std::fmod(dist_v[1], 1);
//            dist_v[2] = std::fmod(dist_v[2], 1);
//            if(std::abs(dist_v[0]) < std::numeric_limits<float>::epsilon()){
//                crit_v[0] = 0;
//            }else{
//                crit_v[0] = dist_v[0] < 0 ? -1 : 1;
//            }
//            if(std::abs(dist_v[1]) < std::numeric_limits<float>::epsilon()){
//                crit_v[1] = 0;
//            }else{
//                crit_v[1] = dist_v[1] < 0 ? -1 : 1;
//            }
//            if(std::abs(dist_v[2]) < std::numeric_limits<float>::epsilon()){
//                crit_v[2] = 0;
//            }else{
//                crit_v[2] = dist_v[2] < 0 ? -1 : 1;
//            }
//            if(!(crit_v[0]||crit_v[1]||crit_v[2])){
//                // TODO: fail here? set flag? overlapping atoms!
//                continue;
//            }
//            dist_v = dist_v * cellvec * celldim;
//            // 0-vector
//            checkBond(i, j, effcut, dist_v, diff_v);
//            if(crit_v[0]){
//                // x, -x
//                checkBond(i, j, effcut, dist_v-crit_v[0]*x,
//                          {{diff_v[0]+crit_v[0],diff_v[1],diff_v[2]}});
//            }
//            if(crit_v[1]){
//                // y, -y
//                checkBond(i, j, effcut, dist_v-crit_v[1]*y,
//                          {{diff_v[0],diff_v[1]+crit_v[1],diff_v[2]}});
//                if(crit_v[0]){
//                    if(crit_v[0] == crit_v[1]){
//                        // x+y, -x-y
//                        checkBond(i, j, effcut, dist_v-crit_v[0]*xy,
//                                  {{diff_v[0]+crit_v[0],
//                                    diff_v[1]+crit_v[1],
//                                    diff_v[2]}});
//                    }else{
//                        // x-y, -x+y
//                        checkBond(i, j, effcut, dist_v-crit_v[0]*xmy,
//                                  {{diff_v[0]+crit_v[0],
//                                    diff_v[1]+crit_v[1],
//                                    diff_v[2]}});
//                    }
//                }
//            }
//            if(crit_v[2]){
//                // z, -z
//                checkBond(i, j, effcut, dist_v-crit_v[2]*z,
//                          {{diff_v[0],diff_v[1],diff_v[2]+crit_v[2]}});
//                if(crit_v[0]){
//                    if(crit_v[0] == crit_v[2]){
//                        // x+z, -x-z
//                        checkBond(i, j, effcut, dist_v-crit_v[0]*xz,
//                                  {{diff_v[0]+crit_v[0],
//                                    diff_v[1],
//                                    diff_v[2]+crit_v[2]}});
//                    }else{
//                        // x-z, -x+z
//                        checkBond(i, j, effcut, dist_v-crit_v[0]*xmz,
//                                  {{diff_v[0]+crit_v[0],
//                                    diff_v[1],
//                                    diff_v[2]+crit_v[2]}});
//                    }
//                }
//                if(crit_v[1]){
//                    if(crit_v[1] == crit_v[2]){
//                        // y+z, -y-z
//                        checkBond(i, j, effcut, dist_v-crit_v[1]*yz,
//                                  {{diff_v[0],
//                                    diff_v[1]+crit_v[1],
//                                    diff_v[2]+crit_v[2]}});
//                    }else{
//                        // y-z, -y+z
//                        checkBond(i, j, effcut, dist_v-crit_v[1]*ymz,
//                                  {{diff_v[0],
//                                    diff_v[1]+crit_v[1],
//                                    diff_v[2]+crit_v[2]}});
//                    }
//                    if(crit_v[0]){
//                        if(crit_v[0] == crit_v[1]){
//                            if(crit_v[0] == crit_v[2]){
//                                // x+y+z, -x-y-z
//                                checkBond(i, j, effcut, dist_v-crit_v[0]*xyz,
//                                          {{diff_v[0]+crit_v[0],
//                                            diff_v[1]+crit_v[1],
//                                            diff_v[2]+crit_v[2]}});
//                            }else{
//                                // x+y-z, -x-y+z
//                                checkBond(i, j, effcut, dist_v-crit_v[0]*xymz,
//                                          {{diff_v[0]+crit_v[0],
//                                            diff_v[1]+crit_v[1],
//                                            diff_v[2]+crit_v[2]}});
//                            }
//                        }else{
//                            if(crit_v[0] == crit_v[2]){
//                                // x-y+z, -x+y-z
//                                checkBond(i, j, effcut, dist_v-crit_v[0]*xmyz,
//                                          {{diff_v[0]+crit_v[0],
//                                            diff_v[1]+crit_v[1],
//                                            diff_v[2]+crit_v[2]}});
//                            }else{
//                                // x-y-z, -x+y+z
//                                checkBond(i, j, effcut, dist_v-crit_v[1]*mxyz,
//                                          {{diff_v[0]+crit_v[0],
//                                            diff_v[1]+crit_v[1],
//                                            diff_v[2]+crit_v[2]}});
//                            }
//                        }
//                    }
//                }
//            }
////            Vec ctemp;
////            ctemp[0] = dist_v[0]*invvec[0][0] + dist_v[1]*invvec[1][0] + dist_v[2]*invvec[2][0];
////            ctemp[1] = dist_v[0]*invvec[0][1] + dist_v[1]*invvec[1][1] + dist_v[2]*invvec[2][1];
////            ctemp[2] = dist_v[0]*invvec[0][2] + dist_v[1]*invvec[1][2] + dist_v[2]*invvec[2][2];
////            dist_v.swap(ctemp);
////            dist_v /= celldim;
////            dist_v += 0.5;
////            int xdiff = std::floor(dist_v[0]);
////            int ydiff = std::floor(dist_v[1]);
////            int zdiff = std::floor(dist_v[2]);
////            dist_v[0] = std::abs(std::fmod(dist_v[0], 1));
////            dist_v[1] = std::abs(std::fmod(dist_v[1], 1));
////            dist_v[2] = std::abs(std::fmod(dist_v[2], 1));
////            dist_v -= 0.5;
////            ctemp[0] = dist_v[0]*cellvec[0][0] + dist_v[1]*cellvec[1][0] + dist_v[2]*cellvec[2][0];
////            ctemp[1] = dist_v[0]*cellvec[0][1] + dist_v[1]*cellvec[1][1] + dist_v[2]*cellvec[2][1];
////            ctemp[2] = dist_v[0]*cellvec[0][2] + dist_v[1]*cellvec[1][2] + dist_v[2]*cellvec[2][2];
////            dist_v.swap(ctemp);
////            dist_v *= celldim;
////            float dist_n = Vec_length(dist_v);
////            if ((0.57 < dist_n) && (dist_n < effcut)){
////                bonds.push_back({i, j, std::sqrt(dist_n), xdiff, ydiff, zdiff});
////            }
//        }
//    }
//    bonds_outdated = false;
//    bondcut_factor = cutfac;
//    bonds_level = BondLevel::Cell;
//}
