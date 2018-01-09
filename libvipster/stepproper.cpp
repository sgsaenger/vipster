#include "stepproper.h"

using namespace Vipster;

//TODO: think about synchronization points for formatted caches

constexpr StepFormatter StepProper::* StepProper::fmtmap[];

StepProper::StepProper(std::shared_ptr<PseMap> pse, AtomFmt fmt, std::string comment)
    : Step{pse, fmt},
      asAlat{this, AtomFmt::Alat},
      asAngstrom{this, AtomFmt::Angstrom},
      asBohr{this, AtomFmt::Bohr},
      asCrystal{this, AtomFmt::Crystal},
      comment{comment}
{
    at_coord = (this->*(fmtmap[(int)fmt])).at_coord;
}

//TODO: move-assignment, copy-and-swap-idiom, constructable from Step&

StepProper::StepProper(const StepProper &rhs)
    : Step{rhs.pse, rhs.fmt},
      asAlat{this, rhs.asAlat},
      asAngstrom{this, rhs.asAngstrom},
      asBohr{this, rhs.asBohr},
      asCrystal{this, rhs.asCrystal},
      at_name{rhs.at_name}, at_charge{rhs.at_charge},
      at_fix{rhs.at_fix}, at_hidden{rhs.at_hidden},
      celldimB{rhs.celldimB}, celldimA{rhs.celldimA},
      cellvec{rhs.cellvec}, invvec{rhs.invvec},
      comment{rhs.comment}
{
    at_coord = (this->*(fmtmap[(int)fmt])).at_coord;
}

StepProper& StepProper::operator=(const StepProper& rhs)
{
    pse = rhs.pse;
    fmt = rhs.fmt;
    comment = rhs.comment;
    asAlat = StepFormatter{this, rhs.asAlat};
    asBohr = StepFormatter{this, rhs.asBohr};
    asAngstrom = StepFormatter{this, rhs.asAngstrom};
    asCrystal = StepFormatter{this, rhs.asCrystal};
    setFmt(fmt);
    return *this;
}

void StepProper::evaluateChanges() const
{
    auto cft = fmtmap[(int)fmt];
    // if there are (only) local changes, invalidate format-caches
    if(at_changed){
        for(auto sft: fmtmap) {
            if ((this->*sft).at_changed) {
                throw Error("Concurrent modification of a "
                            "single Step not supported.");
            } else if (sft == cft) {
                (this->*sft).at_outdated = false;
            } else {
                (this->*sft).at_outdated = true;
            }
        }
        at_changed = false;
        //TODO: mark bonds outdated and similar stuff.
    } else {
        // if there are changes in (only one) formatter,
        // pull them in and invalidate the others
        StepFormatter StepProper::* nft = nullptr;
        for(auto sft: fmtmap) {
            if ((this->*sft).at_changed) {
                if (nft == nullptr) {
                    nft = sft;
                } else {
                    throw Error("Concurrent modification of a "
                                "single Step not supported.");
                }
            }
        }
        if(nft != nullptr) {
            *at_coord = formatAll(*((this->*nft).at_coord), (this->*nft).fmt, fmt);
            (this->*nft).at_changed = false;
            (this->*nft).at_outdated = false;
            for(auto sft: fmtmap) {
                if (sft != nft && sft != cft) {
                    (this->*sft).at_outdated = true;
                }
            }
        }
    }
}

void StepProper::setFmt(AtomFmt format, bool scale)
{
    evaluateChanges();
    if(scale){
        auto oft = fmtmap[(int)fmt];
        auto nft = fmtmap[(int)format];
        for (auto sft: fmtmap){
            if (sft == oft) {
                // old formatter needs new data
                (this->*sft).at_coord = std::make_shared<std::vector<Vec>>();
                (this->*sft).at_outdated = true;
            } else if (sft == nft) {
                // new formatter points to data of old formatter
                (this->*sft).at_coord = at_coord;
                (this->*sft).at_outdated = false;
            } else {
                // other formatters will be renewed when needed
                (this->*sft).at_outdated = true;
            }
        }
    }else{
        auto nft = fmtmap[(int)format];
        // make sure that formatted data is up to date
        (this->*nft).evaluateChanges();
        // let local coordinates point to new formatter
        at_coord = (this->*nft).at_coord;
    }
    fmt = format;
}

StepFormatter& StepProper::asFmt(AtomFmt fmt)
{
    return (this->*StepProper::fmtmap[(int)fmt]);
}

void StepProper::newAtom()
{
    newAtom(AtomProper{});
}

void StepProper::newAtom(const Atom &at)
{
    evaluateChanges();
    at_name.push_back(at.name);
    at_coord->push_back(at.coord);
    at_charge.push_back(at.charge);
    at_fix.push_back(at.fix);
    at_hidden.push_back(at.hidden);
    at_changed = true;
}

void StepProper::newAtoms(size_t i)
{
    evaluateChanges();
    size_t oldNat = getNat();
    size_t nat = oldNat + i;
    at_name.resize(nat);
    at_coord->resize(nat);
    at_charge.resize(nat);
    at_fix.resize(nat);
    at_hidden.resize(nat);
    at_changed = true;
}

void StepProper::delAtom(size_t idx)
{
    evaluateChanges();
    at_name.erase(at_name.begin()+idx);
    at_coord->erase(at_coord->begin()+idx);
    at_charge.erase(at_charge.begin()+idx);
    at_fix.erase(at_fix.begin()+idx);
    at_hidden.erase(at_hidden.begin()+idx);
    at_changed = true;
}

Atom StepProper::operator[](size_t idx)
{
    evaluateChanges();
    return {&at_name[idx], &(*at_coord)[idx], &at_charge[idx],
            &at_fix[idx], &at_hidden[idx], &at_changed};
}

const Atom StepProper::operator[](size_t idx) const
{
    evaluateChanges();
    return {&at_name[idx], &(*at_coord)[idx], &at_charge[idx],
            &at_fix[idx], &at_hidden[idx], &at_changed};
}

void StepProper::setCellDim(float cdm, CdmFmt fmt, bool scale)
{
    if(!(cdm>0))throw Error("Step::setCellDim(): "
                            "cell-dimension needs to be positive");
    evaluateChanges();
    if(scale)
    {
        float ratio = cdm / getCellDim(fmt);
        for(auto& c:*at_coord) {c *= ratio;}
    }
    switch(fmt){
    case CdmFmt::Bohr:
        celldimB = cdm;
        celldimA = cdm*bohrrad;
        break;
    case CdmFmt::Angstrom:
        celldimA = cdm;
        celldimB = cdm*invbohr;
        break;
    }
    at_changed = true;
}

float StepProper::getCellDim(CdmFmt fmt) const noexcept
{
    if (fmt == CdmFmt::Bohr) {
        return celldimB;
    } else {
        return celldimA;
    }
}

void StepProper::setCellVec(const Mat &mat, bool scale)
{
    Mat inv = Mat_inv(mat);
    if (scale) {
        // crystal stays the same
        // other will be updated from updated crystal
        asCrystal.evaluateChanges();
        cellvec = mat;
        invvec = inv;
        if (fmt != AtomFmt::Crystal) {
            *at_coord = formatAll(*(asCrystal.at_coord), AtomFmt::Crystal, fmt);
        }
        at_changed = true;
    } else {
        // keep non-crystal same, invalidate crystal
        evaluateChanges();
        if (fmt != AtomFmt::Crystal) {
            asCrystal.at_outdated = true;
            cellvec = mat;
            invvec = inv;
        } else {
            bool updated{false};
            // if one cache is valid, use it
            for(auto nft:fmtmap) {
                if (nft == &StepProper::asCrystal) continue;
                if ((this->*nft).at_outdated == false) {
                    cellvec = mat;
                    invvec = inv;
                    *at_coord = formatAll(*((this->*nft).at_coord),
                                          (this->*nft).fmt, AtomFmt::Crystal);
                    updated = true;
                    break;
                }
            }
            // else use Alat
            if (!updated) {
                asAlat.evaluateChanges();
                cellvec = mat;
                invvec = inv;
                *at_coord = formatAll(*(asAlat.at_coord),
                                      asAlat.fmt, AtomFmt::Crystal);
            }
            at_changed = true;
        }
    }
}

Mat StepProper::getCellVec() const noexcept
{
    return cellvec;
}

Mat StepProper::getInvVec() const noexcept
{
    return invvec;
}

Vec StepProper::getCenter(CdmFmt fmt, bool com) const noexcept
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
        return formatVec((min+max)/2, this->fmt, (AtomFmt)fmt);
    }else{
        return (cellvec[0]+cellvec[1]+cellvec[2]) * getCellDim(fmt) / 2;
    }
}

void StepProper::setComment(const std::string &s)
{
    comment = s;
}

const std::string& StepProper::getComment() const noexcept
{
    return comment;
}

//const std::vector<Bond>& StepInterface::getBonds() const
//{
//    return getBonds(*bondcut_factor);
//}

//const std::vector<Bond>& StepInterface::getBonds(float cutfac) const
//{
//    //TODO
////    if(bonds_outdated or (cutfac!=bondcut_factor) or (bonds_level<BondLevel::Molecule))
////    {
////        Step::setBonds(cutfac);
////    }
//    return *bonds;
//}

//const std::vector<Bond>& StepInterface::getBondsCell() const
//{
//    return getBondsCell(*bondcut_factor);
//}

//const std::vector<Bond>& StepInterface::getBondsCell(float cutfac) const
//{
//    //TODO
////    if(bonds_outdated or (cutfac!=bondcut_factor) or (bonds_level<BondLevel::Cell))
////    {
////        Step::setBondsCell(cutfac);
////    }
//    return *bonds;
//}

//size_t StepInterface::getNbond() const noexcept
//{
//    return bonds->size();
//}

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

//void StepInterface::checkBond(std::size_t i, std::size_t j, float cutfac, Vec dist, std::array<int, 3> offset) const
//{
//    if (dist[0] > cutfac) return;
//    if (dist[1] > cutfac) return;
//    if (dist[2] > cutfac) return;
//    float dist_n = Vec_dot(dist, dist);
//    if((0.57 < dist_n) && (dist_n < cutfac*cutfac)){
//        bonds->push_back({i, j, std::sqrt(dist_n), offset[0], offset[1], offset[2]});
//    }
//}

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
