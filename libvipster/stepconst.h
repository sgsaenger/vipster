#ifndef STEPCONST_H
#define STEPCONST_H

#include "atom.h"
#include "bond.h"
#include "cell.h"
#include "vec.h"
#include "settings.h"
#include "stepsel.h"

#include <memory>
#include <vector>
#include <set>
#include <functional>
#include <algorithm>

namespace Vipster {

/*
 * Base for all Step-like containers
 *
 * Implements const-interface for common state:
 * - Atoms (atomcontainer must be provided as template argument)
 * - Bonds
 * - Cell
 * - Comment
 * and provides the storage
 */
template<typename T>
class StepConst
{
public:
    virtual ~StepConst() = default;

    using source = T;
    // 'mutable' iterator must be defined for Selection-template
    using iterator = typename T::const_iterator;
    using const_iterator = typename T::const_iterator;
    using reverse_iterator = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;
    using constSelection = SelectionBase<StepConst, T>;

    void evaluateCache() const; // must be implemented in specialization

    // Don't know how to mask this yet
    std::shared_ptr<PseMap> pse;

    // Selection
    constSelection select(std::string filter) const
    {
        return {pse, at_fmt, this, filter, bonds, cell, comment};
    }
    constSelection select(SelectionFilter filter) const
    {
        return {pse, at_fmt, this, filter, bonds, cell, comment};
    }

    // Atoms
    size_t          getNat() const noexcept; // must be implemented in specialization
    constAtom       operator[](size_t i) const noexcept
    {
        return *const_iterator{atoms, at_fmt, i};
    }
    const source&   getAtoms() const noexcept
    {
        return *atoms;
    }
    const_iterator   begin() const noexcept
    {
        return const_iterator{atoms, at_fmt, 0};
    }
    const_iterator   cbegin() const noexcept
    {
        return const_iterator{atoms, at_fmt, 0};
    }
    const_iterator   end() const noexcept
    {
        return const_iterator{atoms, at_fmt, getNat()};
    }
    const_iterator   cend() const noexcept
    {
        return const_iterator{atoms, at_fmt, getNat()};
    }
    const_reverse_iterator rbegin() const noexcept
    {
        return std::make_reverse_iterator(end());
    }
    const_reverse_iterator crbegin() const noexcept
    {
        return std::make_reverse_iterator(cend());
    }
    const_reverse_iterator rend() const noexcept
    {
        return std::make_reverse_iterator(begin());
    }
    const_reverse_iterator crend() const noexcept
    {
        return std::make_reverse_iterator(cbegin());
    }

    // Comment
    const std::string&  getComment() const noexcept
    {
        return *comment;
    }

    // Types
    std::set<std::string>   getTypes() const
    {
        std::set<std::string> tmp{};
        for(const auto& at: *this){
            tmp.insert(at.name);
        }
        return tmp;
    }
    size_t                  getNtyp() const
    {
        return getTypes().size();
    }

    // Format
    AtomFmt             getFmt() const noexcept
    {
        return at_fmt;
    }
    void                setFmt(AtomFmt tgt) const
    {
        if(tgt == at_fmt){ return; }
        asFmt(tgt).evaluateCache();
        at_fmt = tgt;
    }
    StepConst           asFmt(AtomFmt tgt) const
    {
        auto tmp = StepConst{pse, tgt, atoms, bonds, cell, comment};
        tmp.evaluateCache();
        return tmp;
    }
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
    Vec                 formatVec(Vec in, AtomFmt source, AtomFmt target) const
    {
        return getFormatter(source, target)(in);
    }
    std::vector<Vec>    formatAll(std::vector<Vec> in, AtomFmt source,
                                  AtomFmt target) const
    {
        if ((source == target) || in.empty()){
            return in;
        }
        auto op = getFormatter(source, target);
        std::transform(in.begin(), in.end(), in.begin(), op);
        return in;
    }

    // Bonds
    const std::vector<Bond>&    getBonds(float cutfac=settings.bondCutFac.val,
                                         BondLevel l=settings.bondLvl.val,
                                         BondFrequency update=settings.bondFreq.val) const
    {
        if(cutfac < 0){
            throw Error("bond cutoff factor must be positive");
        }
        if(((update == BondFrequency::Always) or
            ((update == BondFrequency::Once) and not bonds->setOnce))
           and
           (bonds->outdated or
            (!float_comp(cutfac, bonds->cutoff_factor)) or
            (bonds->level != l)))
        {
            setBonds(l, cutfac);
        }
        return bonds->bonds;
    }
    size_t                      getNbond() const
    {
        return getBonds().size();
    }
    void                        setBonds(BondLevel l=settings.bondLvl.val,
                                         float cutfac=settings.bondCutFac.val) const
    {
        if(cutfac < 0){
            throw Error("bond cutoff factor must be positive");
        }
        bonds->bonds.clear();
        if(!cutfac){
            return;
        }
        if(!getNat()){
            l = BondLevel::None;
        }
        if((l==BondLevel::Cell) && !cell->enabled){ l = BondLevel::Molecule; }
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
        bonds->setOnce = true;
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
        }
        return cell->dimAngstrom;
    }
    void    enableCell(bool val) const noexcept
    {
        cell->enabled = val;
    }
    Mat     getCellVec() const noexcept
    {
        return cell->cellvec;
    }
    Vec     getCom() const noexcept
    {
        return getCom(at_fmt);
    }
    Vec     getCom(AtomFmt fmt) const noexcept
    {
        if(!getNat()){
            return Vec{{0,0,0}};
        }
        Vec min{{std::numeric_limits<float>::max(),
                 std::numeric_limits<float>::max(),
                 std::numeric_limits<float>::max()}};
        Vec max{{std::numeric_limits<float>::min(),
                 std::numeric_limits<float>::min(),
                 std::numeric_limits<float>::min()}};
        for(const auto& at: *this){
            min[0]=std::min(min[0],at.coord[0]);
            min[1]=std::min(min[1],at.coord[1]);
            min[2]=std::min(min[2],at.coord[2]);
            max[0]=std::max(max[0],at.coord[0]);
            max[1]=std::max(max[1],at.coord[1]);
            max[2]=std::max(max[2],at.coord[2]);
        }
        return formatVec((min+max)/2, at_fmt, fmt);
    }
    Vec     getCenter(CdmFmt fmt, bool com=false) const noexcept
    {
        if(com || !cell->enabled){
            return getCom(static_cast<AtomFmt>(fmt));
        }
        const Mat& cv = cell->cellvec;
        return (cv[0]+cv[1]+cv[2]) * getCellDim(fmt) / 2;
    }

protected:
    StepConst(std::shared_ptr<PseMap> pse, AtomFmt fmt,
             std::shared_ptr<source> atoms, std::shared_ptr<BondList> bonds,
             std::shared_ptr<CellData> cell, std::shared_ptr<std::string> comment)
        : pse{pse}, at_fmt{fmt}, atoms{atoms}, bonds{bonds}, cell{cell}, comment{comment}
    {}
    // Data
    mutable AtomFmt                 at_fmt; // mutable, only controls visibility
    std::shared_ptr<source>         atoms;  // immutable
    std::shared_ptr<BondList>       bonds;  // mutable, is only cache
    std::shared_ptr<CellData>       cell;   // immutable except enabled/disable-state
    std::shared_ptr<std::string>    comment; // immutable

private:

    // Bonds
    void    setBondsMolecule(float cutfac) const
    {
        const AtomFmt fmt = (this->at_fmt == AtomFmt::Angstrom) ? AtomFmt::Angstrom : AtomFmt::Bohr;
        const float fmtscale{(fmt == AtomFmt::Angstrom) ? invbohr : 1};
        auto tgtFmt = asFmt(fmt);
        tgtFmt.evaluateCache();
        std::vector<Bond>& bonds = this->bonds->bonds;
        auto at_i = tgtFmt.begin();
        for (auto at_i=tgtFmt.begin(); at_i!=tgtFmt.end(); ++at_i)
        {
            float cut_i = at_i->pse->bondcut;
            if (cut_i<0){
                continue;
            }
            for (auto at_j=tgtFmt.begin()+at_i.getIdx()+1; at_j != tgtFmt.end(); ++at_j){
                float cut_j = at_j->pse->bondcut;
                if (cut_j<0) {
                    continue;
                }
                float effcut = (cut_i + cut_j) * cutfac;
                Vec dist_v = at_i->coord - at_j->coord;
                if (((dist_v[0] *= fmtscale) > effcut) ||
                    ((dist_v[1] *= fmtscale) > effcut) ||
                    ((dist_v[2] *= fmtscale) > effcut)) {
                    continue;
                }
                float dist_n = Vec_dot(dist_v, dist_v);
                if((0.57f < dist_n) && (dist_n < effcut*effcut)) {
                    bonds.push_back({at_i.getIdx(), at_j.getIdx(), std::sqrt(dist_n), 0, 0, 0});
                }
            }
        }
        this->bonds->outdated = false;
        this->bonds->level = BondLevel::Molecule;
    }

    void setBondsCell(float cutfac) const
    {
        auto& bonds = this->bonds->bonds;
        const auto asCrystal = asFmt(AtomFmt::Crystal);
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
        auto checkBond = [&](std::size_t i, std::size_t j, float effcut,
                             const Vec& dist, const std::array<int16_t, 3>& offset)
        {
            if ((dist[0]>effcut) || (dist[1]>effcut) || (dist[2]>effcut)) {
                return;
            }
            float dist_n = Vec_dot(dist, dist);
            if ((0.57f < dist_n) && (dist_n < effcut*effcut)) {
                bonds.push_back({i, j, std::sqrt(dist_n), offset[0], offset[1], offset[2]});
            }
        };
        for(auto at_i = asCrystal.begin(); at_i != asCrystal.end(); ++at_i){
            size_t i = at_i.getIdx();
            float cut_i = at_i->pse->bondcut;
            if (cut_i<0) {
                continue;
            }
            for(auto at_j = at_i+1; at_j != asCrystal.end(); ++at_j){
                size_t j = at_j.getIdx();
                float cut_j = at_j->pse->bondcut;
                if (cut_j<0) {
                    continue;
                }
                float effcut = (cut_i + cut_j) * cutfac;
                Vec dist_v = at_i->coord - at_j->coord;
                // diff_v contains integer distance in cell-units
                std::transform(dist_v.begin(), dist_v.end(), diff_v.begin(), truncf);
                // dist_v now contains distance inside of cell
                std::transform(dist_v.begin(), dist_v.end(), dist_v.begin(),
                    [](float f){return std::fmod(f,1);});
                // crit_v contains direction of distance-vector
                std::transform(dist_v.begin(), dist_v.end(), crit_v.begin(),
                    [](float f){
                        return (std::abs(f) < std::numeric_limits<float>::epsilon())?
                                    0 : ((f<0) ? -1 : 1);
                    });
                if(!((crit_v[0] != 0)||(crit_v[1] != 0)||(crit_v[2] != 0))){
                    // TODO: fail here? set flag? overlapping atoms!
                    continue;
                }
                // convert dist_v to bohr
                dist_v = dist_v * cell->cellvec * cell->dimBohr;
                // 0-vector
                checkBond(i, j, effcut, dist_v, diff_v);
                if(crit_v[0] != 0){
                    // x, -x
                    checkBond(i, j, effcut, dist_v-crit_v[0]*x,
                              {{static_cast<int16_t>(diff_v[0]+crit_v[0]),diff_v[1],diff_v[2]}});
                }
                if(crit_v[1] != 0){
                    // y, -y
                    checkBond(i, j, effcut, dist_v-crit_v[1]*y,
                              {{diff_v[0],static_cast<int16_t>(diff_v[1]+crit_v[1]),diff_v[2]}});
                    if(crit_v[0] != 0){
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
                if(crit_v[2] != 0){
                    // z, -z
                    checkBond(i, j, effcut, dist_v-crit_v[2]*z,
                              {{diff_v[0],diff_v[1],static_cast<int16_t>(diff_v[2]+crit_v[2])}});
                    if(crit_v[0] != 0){
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
                    if(crit_v[1] != 0){
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
                        if(crit_v[0] != 0){
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
            }
        }
        this->bonds->outdated = false;
        this->bonds->level = BondLevel::Cell;
    }
};

}

#endif // STEPCONST_H
