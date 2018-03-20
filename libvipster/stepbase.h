#ifndef BONDCELLINTERFACE_H
#define BONDCELLINTERFACE_H

#include "atom.h"
#include "bond.h"
#include "cell.h"
#include "vec.h"
#include "config.h"

#include <memory>
#include <vector>
#include <set>
#include <functional>
#include <algorithm>

namespace Vipster {

/*
 * Base for all Step-like containers
 *
 * Uses CRTP to access Atoms and cache
 * Implements interface for Bonds, Cell, Comment
 */
template<typename T>
class StepBase
{
public:
    virtual ~StepBase()=default;

    // Don't know how to mask this yet
    std::shared_ptr<PseMap> pse;

    // Comment
    void                setComment(const std::string& s)
    {
        *comment = s;
    }
    const std::string&  getComment() const noexcept
    {
        return *comment;
    }

    // Types
    std::set<std::string>   getTypes() const
    {
        //TODO: use PseEntries instead of string?
        std::set<std::string> set;
        for(auto& at: *static_cast<const T*>(this)){
            set.insert(at.name);
        }
        return set;
    }
    size_t                  getNtyp() const
    {
        return getTypes().size();
    }

    // Format
    AtomFmt getFmt() const noexcept{
        return at_fmt;
    }
    virtual T&          asFmt(AtomFmt)=0;
    virtual const T&    asFmt(AtomFmt) const=0;

    // Bonds
    const std::vector<Bond>&    getBonds(BondLevel l=BondLevel::Cell,
                                         bool update=true) const
    {
        return getBonds(bonds->cutoff_factor, l, update);
    }
    const std::vector<Bond>&    getBonds(float cutfac,
                                         BondLevel l=BondLevel::Cell,
                                         bool update=true) const
    {
        evaluateCache();
        if((l==BondLevel::Cell) && !cell->enabled){ l = BondLevel::Molecule; }
        if(update and (bonds->outdated or
                       (cutfac != bonds->cutoff_factor) or
                       (bonds->level < l)))
        {
            setBonds(l, cutfac);
        }
        return bonds->bonds;
    }
    size_t                      getNbond() const
    {
        return getBonds().size();
    }
    void    setBonds(BondLevel l, float cutfac) const
    {
        bonds->bonds.clear();
        if(!static_cast<const T*>(this)->getNat()) l = BondLevel::None;
        switch(l){
        case BondLevel::None:
            break;
        case BondLevel::Molecule:
            setBondsMolecule(cutfac);
            break;
        case BondLevel::Cell:
            setBondsCell(cutfac);
            break;
        }
        bonds->outdated = false;
        bonds->cutoff_factor = cutfac;
        bonds->level = l;
    }

    // Cell
    bool    hasCell() const noexcept
    {
        return cell->enabled;
    }
    float   getCellDim(CdmFmt fmt) const noexcept
    {
        if (fmt == CdmFmt::Bohr) {
            return cell->dimBohr;
        }else{
            return cell->dimAngstrom;
        }
    }
    Mat     getCellVec() const noexcept
    {
        return cell->cellvec;
    }
    Vec     getCenter(CdmFmt fmt, bool com=false) const noexcept
    {
        if((com && static_cast<const T*>(this)->getNat()) || !cell->enabled){
            Vec min{{std::numeric_limits<float>::max(),
                     std::numeric_limits<float>::max(),
                     std::numeric_limits<float>::max()}};
            Vec max{{std::numeric_limits<float>::min(),
                     std::numeric_limits<float>::min(),
                     std::numeric_limits<float>::min()}};
            for(auto& at:*static_cast<const T*>(this)){
                min[0]=std::min(min[0],at.coord[0]);
                min[1]=std::min(min[1],at.coord[1]);
                min[2]=std::min(min[2],at.coord[2]);
                max[0]=std::max(max[0],at.coord[0]);
                max[1]=std::max(max[1],at.coord[1]);
                max[2]=std::max(max[2],at.coord[2]);
            }
            return formatVec((min+max)/2, at_fmt, static_cast<AtomFmt>(fmt));
        }else if(cell->enabled){
            const Mat& cv = cell->cellvec;
            return (cv[0]+cv[1]+cv[2]) * getCellDim(fmt) / 2;
        }else{
            return Vec{{0,0,0}};
        }
    }
protected:
    virtual void evaluateCache()const =0;
    StepBase(std::shared_ptr<PseMap> pse, AtomFmt fmt, std::shared_ptr<BondList> bonds,
             std::shared_ptr<CellData> cell, std::shared_ptr<std::string> comment)
        : pse{pse}, at_fmt{fmt}, bonds{bonds}, cell{cell}, comment{comment} {}
    // Data
    AtomFmt                         at_fmt;
    std::shared_ptr<BondList>       bonds;
    std::shared_ptr<CellData>       cell;
    std::shared_ptr<std::string>    comment;
    // Format
    std::function<Vec(const Vec&)>  getFormatter(AtomFmt source,
                                                 AtomFmt target) const
    {
        float fac{};
        Mat fmat{};
        switch(source) {
        case AtomFmt::Bohr:
            switch(target){
            case AtomFmt::Angstrom:
                return [](const Vec& v){return v*Vipster::bohrrad;};
            case AtomFmt::Alat:
                fac = 1/getCellDim(CdmFmt::Bohr);
                return [fac](const Vec& v){return v*fac;};
            case AtomFmt::Crystal:
                fmat = Mat_inv(getCellVec())/getCellDim(CdmFmt::Bohr);
                return [fmat](const Vec& v){return v*fmat;};
            default:
                break;
            }
            break;
        case AtomFmt::Angstrom:
            switch(target){
            case AtomFmt::Bohr:
                return [](const Vec& v){return v*Vipster::invbohr;};
            case AtomFmt::Alat:
                fac = 1/getCellDim(CdmFmt::Angstrom);
                return [fac](const Vec& v){return v*fac;};
            case AtomFmt::Crystal:
                fmat = Mat_inv(getCellVec())/getCellDim(CdmFmt::Angstrom);
                return [fmat](const Vec& v){return v*fmat;};
            default:
                break;
            }
            break;
        case AtomFmt::Alat:
            switch(target){
            case AtomFmt::Angstrom:
                fac = getCellDim(CdmFmt::Angstrom);
                return [fac](const Vec& v){return v*fac;};
            case AtomFmt::Bohr:
                fac = getCellDim(CdmFmt::Bohr);
                return [fac](const Vec& v){return v*fac;};
            case AtomFmt::Crystal:
                fmat = Mat_inv(getCellVec());
                return [fmat](const Vec& v){return v*fmat;};
            default:
                break;
            }
            break;
        case AtomFmt::Crystal:
            switch(target){
            case AtomFmt::Angstrom:
                fmat = getCellVec()*getCellDim(CdmFmt::Angstrom);
                return [fmat](const Vec& v){return v*fmat;};
            case AtomFmt::Alat:
                fmat = getCellVec();
                return [fmat](const Vec& v){return v*fmat;};
            case AtomFmt::Bohr:
                fmat = getCellVec()*getCellDim(CdmFmt::Bohr);
                return [fmat](const Vec& v){return v*fmat;};
            default:
                break;
            }
        }
        return [](const Vec& v){return v;};
    }
    Vec                             formatVec(Vec in, AtomFmt source,
                                              AtomFmt target) const
    {
        return getFormatter(source, target)(in);
    }
    std::vector<Vec>                formatAll(std::vector<Vec> in,
                                              AtomFmt source,
                                              AtomFmt target) const
    {
        if ((source == target) || (in.size() == 0)) return in;
        auto op = getFormatter(source, target);
        std::transform(in.begin(), in.end(), in.begin(), op);
        return in;
    }

private:
    // Bonds
    void    setBondsMolecule(float cutfac) const
    {
        AtomFmt fmt = (this->at_fmt == AtomFmt::Angstrom) ? AtomFmt::Angstrom : AtomFmt::Bohr;
        float fmtscale{(fmt == AtomFmt::Angstrom) ? invbohr : 1};
        const T& asFmt = static_cast<const T*>(this)->asFmt(fmt);
        asFmt.evaluateCache();
        std::vector<Bond>& bonds = this->bonds->bonds;
        auto at_i = asFmt.begin();
        for (auto at_i=asFmt.begin(); at_i!=asFmt.end(); ++at_i)
        {
            float cut_i = (*at_i->pse).bondcut;
            if (cut_i<0) continue;
            for (auto at_j=asFmt.begin()+at_i.getIdx()+1; at_j != asFmt.end(); ++at_j){
                float cut_j = (*at_j->pse).bondcut;
                if (cut_j<0) continue;
                float effcut = (cut_i + cut_j) * cutfac;
                Vec dist_v = at_i->coord - at_j->coord;
                if (((dist_v[0] *= fmtscale) > effcut) ||
                    ((dist_v[1] *= fmtscale) > effcut) ||
                    ((dist_v[2] *= fmtscale) > effcut)) continue;
                float dist_n = Vec_dot(dist_v, dist_v);
                if((0.57f < dist_n) && (dist_n < effcut*effcut)) {
                    bonds.push_back({at_i.getIdx(), at_j.getIdx(), std::sqrt(dist_n), 0, 0, 0});
                }
            }
        }
        this->bonds->outdated = false;
        this->bonds->level = BondLevel::Molecule;
    }
    void    checkBond(std::size_t i, std::size_t j, float effcut,
                      const Vec& dist, const std::array<int16_t, 3>& offset) const
    {
        auto& bonds = this->bonds->bonds;
        if ((dist[0]>effcut) || (dist[1]>effcut) || (dist[2]>effcut)) return;
        float dist_n = Vec_dot(dist, dist);
        if ((0.57f < dist_n) && (dist_n < effcut*effcut)) {
            bonds.push_back({i, j, std::sqrt(dist_n), offset[0], offset[1], offset[2]});
        }
    }
    void    setBondsCell(float cutfac) const
    {
        const T& asCrystal = asFmt(AtomFmt::Crystal);
        asCrystal.evaluateCache();
        const Vec x = getCellVec()[0] * getCellDim(CdmFmt::Bohr);
        const Vec y = getCellVec()[1] * getCellDim(CdmFmt::Bohr);
        const Vec z = getCellVec()[2] * getCellDim(CdmFmt::Bohr);
        const Vec xy   = x+y;
        const Vec xmy  = x-y;
        const Vec xz   = x+z;
        const Vec xmz  = x-z;
        const Vec yz   = y+z;
        const Vec ymz  = y-z;
        const Vec xyz  = xy + z;
        const Vec xymz = xy - z;
        const Vec xmyz = xz - y;
        const Vec mxyz = yz - x;
        std::array<int16_t, 3> diff_v, crit_v;
        size_t nat = static_cast<const T*>(this)->getNat();
        auto at_i = asCrystal.begin();
        for (size_t i=0; i<nat; ++i) {
            float cut_i = (*at_i->pse).bondcut;
            if (cut_i<0) continue;
            auto at_j = asCrystal.begin();
            for (size_t j=0; j<nat; ++j) {
                float cut_j = (*at_j->pse).bondcut;
                if (cut_j<0) continue;
                float effcut = (cut_i + cut_j) * cutfac;
                Vec dist_v = at_i->coord - at_j->coord;
                std::transform(dist_v.begin(), dist_v.end(), diff_v.begin(), truncf);
                std::transform(dist_v.begin(), dist_v.end(), dist_v.begin(),
                    [](float f){return std::fmod(f,1);});
                std::transform(dist_v.begin(), dist_v.end(), crit_v.begin(),
                    [](float f){
                        return (std::abs(f) < std::numeric_limits<float>::epsilon())?
                                    0 : ((f<0) ? -1 : 1);
                    });
                if(!(crit_v[0]||crit_v[1]||crit_v[2])){
                    // TODO: fail here? set flag? overlapping atoms!
                    continue;
                }
                dist_v = dist_v * cell->cellvec * cell->dimBohr;
                // 0-vector
                checkBond(i, j, effcut, dist_v, diff_v);
                if(crit_v[0]){
                    // x, -x
                    checkBond(i, j, effcut, dist_v-crit_v[0]*x,
                              {{static_cast<int16_t>(diff_v[0]+crit_v[0]),diff_v[1],diff_v[2]}});
                }
                if(crit_v[1]){
                    // y, -y
                    checkBond(i, j, effcut, dist_v-crit_v[1]*y,
                              {{diff_v[0],static_cast<int16_t>(diff_v[1]+crit_v[1]),diff_v[2]}});
                    if(crit_v[0]){
                        if(crit_v[0] == crit_v[1]){
                            // x+y, -x-y
                            checkBond(i, j, effcut, dist_v-crit_v[0]*xy,
                                      {{static_cast<int16_t>(diff_v[0]+crit_v[0]),
                                        static_cast<int16_t>(diff_v[1]+crit_v[1]),
                                        diff_v[2]}});
                        }else{
                            // x-y, -x+y
                            checkBond(i, j, effcut, dist_v-crit_v[0]*xmy,
                                      {{static_cast<int16_t>(diff_v[0]+crit_v[0]),
                                        static_cast<int16_t>(diff_v[1]+crit_v[1]),
                                        diff_v[2]}});
                        }
                    }
                }
                if(crit_v[2]){
                    // z, -z
                    checkBond(i, j, effcut, dist_v-crit_v[2]*z,
                              {{diff_v[0],diff_v[1],static_cast<int16_t>(diff_v[2]+crit_v[2])}});
                    if(crit_v[0]){
                        if(crit_v[0] == crit_v[2]){
                            // x+z, -x-z
                            checkBond(i, j, effcut, dist_v-crit_v[0]*xz,
                                      {{static_cast<int16_t>(diff_v[0]+crit_v[0]),
                                        diff_v[1],
                                        static_cast<int16_t>(diff_v[2]+crit_v[2])}});
                        }else{
                            // x-z, -x+z
                            checkBond(i, j, effcut, dist_v-crit_v[0]*xmz,
                                      {{static_cast<int16_t>(diff_v[0]+crit_v[0]),
                                        diff_v[1],
                                        static_cast<int16_t>(diff_v[2]+crit_v[2])}});
                        }
                    }
                    if(crit_v[1]){
                        if(crit_v[1] == crit_v[2]){
                            // y+z, -y-z
                            checkBond(i, j, effcut, dist_v-crit_v[1]*yz,
                                      {{diff_v[0],
                                        static_cast<int16_t>(diff_v[1]+crit_v[1]),
                                        static_cast<int16_t>(diff_v[2]+crit_v[2])}});
                        }else{
                            // y-z, -y+z
                            checkBond(i, j, effcut, dist_v-crit_v[1]*ymz,
                                      {{diff_v[0],
                                        static_cast<int16_t>(diff_v[1]+crit_v[1]),
                                        static_cast<int16_t>(diff_v[2]+crit_v[2])}});
                        }
                        if(crit_v[0]){
                            if(crit_v[0] == crit_v[1]){
                                if(crit_v[0] == crit_v[2]){
                                    // x+y+z, -x-y-z
                                    checkBond(i, j, effcut, dist_v-crit_v[0]*xyz,
                                              {{static_cast<int16_t>(diff_v[0]+crit_v[0]),
                                                static_cast<int16_t>(diff_v[1]+crit_v[1]),
                                                static_cast<int16_t>(diff_v[2]+crit_v[2])}});
                                }else{
                                    // x+y-z, -x-y+z
                                    checkBond(i, j, effcut, dist_v-crit_v[0]*xymz,
                                              {{static_cast<int16_t>(diff_v[0]+crit_v[0]),
                                                static_cast<int16_t>(diff_v[1]+crit_v[1]),
                                                static_cast<int16_t>(diff_v[2]+crit_v[2])}});
                                }
                            }else{
                                if(crit_v[0] == crit_v[2]){
                                    // x-y+z, -x+y-z
                                    checkBond(i, j, effcut, dist_v-crit_v[0]*xmyz,
                                              {{static_cast<int16_t>(diff_v[0]+crit_v[0]),
                                                static_cast<int16_t>(diff_v[1]+crit_v[1]),
                                                static_cast<int16_t>(diff_v[2]+crit_v[2])}});
                                }else{
                                    // x-y-z, -x+y+z
                                    checkBond(i, j, effcut, dist_v-crit_v[1]*mxyz,
                                              {{static_cast<int16_t>(diff_v[0]+crit_v[0]),
                                                static_cast<int16_t>(diff_v[1]+crit_v[1]),
                                                static_cast<int16_t>(diff_v[2]+crit_v[2])}});
                                }
                            }
                        }
                    }
                }
                ++at_j;
            }
            ++at_i;
        }
        bonds->outdated = false;
        bonds->level = BondLevel::Cell;
    }
};

}

#endif // BONDCELLINTERFACE_H
