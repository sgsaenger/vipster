#ifndef STEPCONST_H
#define STEPCONST_H

#include "atom.h"
#include "bond.h"
#include "cell.h"
#include "vec.h"
#include "data.h"
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
    friend T;
public:
    virtual ~StepConst() = default;

    using source = T;
    // 'mutable' iterator must be defined for Selection-template
    using iterator = typename T::const_iterator;
    using const_iterator = typename T::const_iterator;
    using reverse_iterator = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;
    using constSelection = SelectionBase<StepConst, T>;

    void evaluateCache() const
    {
        atoms->evaluateCache(*this);
    }

    // Don't know how to mask this yet
    std::shared_ptr<PseMap> pse;

    // Selection
    constSelection select(std::string filter) const
    {
        return {pse, at_fmt, this, filter, cell, comment};
    }
    constSelection select(SelectionFilter filter) const
    {
        return {pse, at_fmt, this, filter, cell, comment};
    }

    // Atoms
    size_t          getNat() const noexcept
    {
        return atoms->getNat();
    }
    constAtom       operator[](size_t i) const noexcept
    {
        return *const_iterator{atoms, at_fmt, i};
    }
    constAtom       at(size_t i) const
    {
        if(i>=getNat()){
            throw Error("Atom-index out of bounds");
        }else{
            return operator[](i);
        }
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
        asFmt(tgt);
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
        if((l==BondLevel::Cell) && !cell->enabled){ l = BondLevel::Molecule; }
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
    void setBondsMolecule(float cutfac) const
    {
        if(getNat() >= 1000){
            setBondsMoleculeSplit(cutfac);
        }else{
            setBondsMoleculeTrivial(cutfac);
        }
    }
    void setBondsMoleculeSplit(float cutfac) const
    {
        // get suitable absolute representation
        const AtomFmt fmt = (this->at_fmt == AtomFmt::Angstrom) ? AtomFmt::Angstrom : AtomFmt::Bohr;
        const float fmtscale{(fmt == AtomFmt::Angstrom) ? invbohr : 1};
        auto tgtFmt = asFmt(fmt);
        tgtFmt.evaluateCache();
        // get bounds of system and largest cutoff
        Vec min{0,0,0};
        Vec max{0,0,0};
        float cut{0};
        for (const auto& at:tgtFmt) {
            min[0] = std::min(at.coord[0], min[0]);
            min[1] = std::min(at.coord[1], min[1]);
            min[2] = std::min(at.coord[2], min[2]);
            max[0] = std::max(at.coord[0], max[0]);
            max[1] = std::max(at.coord[1], max[1]);
            max[2] = std::max(at.coord[2], max[2]);
            cut = std::max(at.pse->bondcut, cut);
        }
        // fragment spanned space
        Vec diff = max - min;
        cut = 5*cut;
        Vec size_split = diff;
        SizeVec n_split{1,1,1};
        if(diff[0] >= cut){
            n_split[0] = static_cast<size_t>(std::round(diff[0] / cut));
            size_split[0] = diff[0] / n_split[0];
        }
        if(diff[1] >= cut){
            n_split[1] = static_cast<size_t>(std::round(diff[1] / cut));
            size_split[1] = diff[1] / n_split[1];
        }
        if(diff[2] >= cut){
            n_split[2] = static_cast<size_t>(std::round(diff[2] / cut));
            size_split[2] = diff[2] / n_split[2];
        }
        // put atoms in boxes
        DataGrid3D<std::vector<size_t>> boxes{n_split};
        for(auto it=tgtFmt.begin(); it!=tgtFmt.end(); ++it){
            auto findBox = [&](size_t dir){
                return std::max(size_t{0},
                                std::min(n_split[dir]-1,
                                         static_cast<size_t>((it->coord[dir]-min[dir])/size_split[dir])));
            };
            boxes(findBox(0), findBox(1), findBox(2)).push_back(it.getIdx());
        }
        // get bonds by iterating over boxes and their neighbors
        auto& bonds = this->bonds->bonds;
        for(size_t x=0; x<boxes.extent[0]; ++x){
            for(size_t y=0; y<boxes.extent[1]; ++y){
                for(size_t z=0; z<boxes.extent[2]; ++z){
                    for (auto it_i=boxes(x,y,z).begin(); it_i != boxes(x,y,z).end(); ++it_i) {
                        auto i = *it_i;
                        const auto& at_i = tgtFmt[i];
                        auto cut_i = at_i.pse->bondcut;
                        if (cut_i <= 0){
                            continue;
                        }
                        // loop over rest of current box
                        for (auto it_j = it_i+1; it_j!=boxes(x,y,z).end(); ++it_j) {
                            auto j = *it_j;
                            const auto& at_j = tgtFmt[j];
                            auto cut_j = at_j.pse->bondcut;
                            if (cut_j <= 0) {
                                continue;
                            }
                            auto effcut = (cut_i + cut_j) * cutfac;
                            Vec dist_v = at_i.coord - at_j.coord;
                            if (((dist_v[0] *= fmtscale) > effcut) ||
                                ((dist_v[1] *= fmtscale) > effcut) ||
                                ((dist_v[2] *= fmtscale) > effcut)) {
                                continue;
                            }
                            auto dist_n = Vec_dot(dist_v, dist_v);
                            if((0.57f < dist_n) && (dist_n < effcut*effcut)) {
                                bonds.push_back({i, j, std::sqrt(dist_n), 0, 0, 0});
                            }
                        }
                        // visit neighboring boxes
                        auto bondcheck = [&](size_t xt, size_t yt, size_t zt){
                            for (const auto& j: boxes(xt, yt, zt)) {
                                const auto& at_j = tgtFmt[j];
                                auto cut_j = at_j.pse->bondcut;
                                if (cut_j <= 0) {
                                    continue;
                                }
                                auto effcut = (cut_i + cut_j) * cutfac;
                                Vec dist_v = at_i.coord - at_j.coord;
                                if (((dist_v[0] *= fmtscale) > effcut) ||
                                    ((dist_v[1] *= fmtscale) > effcut) ||
                                    ((dist_v[2] *= fmtscale) > effcut)) {
                                    continue;
                                }
                                auto dist_n = Vec_dot(dist_v, dist_v);
                                if((0.57f < dist_n) && (dist_n < effcut*effcut)) {
                                    bonds.push_back({i, j, std::sqrt(dist_n), 0, 0, 0});
                                }
                            }
                        };
                        // +x (includes -x)
                        if (x < boxes.extent[0]-1) {
                            bondcheck(x+1, y, z);
                            // +x+y (includes -x-y)
                            if (y < boxes.extent[1]-1) {
                                bondcheck(x+1, y+1, z);
                                // +x+y+z (includes -x-y-z)
                                if (z < boxes.extent[2]-1) {
                                    bondcheck(x+1, y+1, z+1);
                                }
                                // +x+y-z (includes -x-y+z)
                                if (z > 0) {
                                    bondcheck(x+1, y+1, z-1);
                                }
                            }
                            // +x-y (includes -x+y)
                            if (y > 0) {
                                bondcheck(x+1, y-1, z);
                                // +x-y+z (includes -x+y-z)
                                if (z < boxes.extent[2]-1) {
                                    bondcheck(x+1, y-1, z+1);
                                }
                                // +x-y-z (includes -x+y+z)
                                if (z > 0) {
                                    bondcheck(x+1, y-1, z-1);
                                }
                            }
                            // +x+z (includes -x-z)
                            if (z < boxes.extent[2]-1) {
                                bondcheck(x+1, y, z+1);
                            }
                            // +x-z (includex -x+z)
                            if (z > 0) {
                                bondcheck(x+1, y, z-1);
                            }
                        }
                        // +y (includes -y)
                        if (y < boxes.extent[1]-1) {
                            bondcheck(x, y+1, z);
                            // +y+z (includes -y-z)
                            if (z < boxes.extent[2]-1) {
                                bondcheck(x, y+1, z+1);
                            }
                            // +y-z (includes -y+z)
                            if (z > 0) {
                                bondcheck(x, y+1, z-1);
                            }
                        }
                        // +z (includes -z)
                        if (z < boxes.extent[2]-1) {
                            bondcheck(x, y, z+1);
                        }
                    }
                }
            }
        }
    }
    void setBondsMoleculeTrivial(float cutfac) const
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
            if (cut_i<=0){
                continue;
            }
            for (auto at_j=tgtFmt.begin()+at_i.getIdx()+1; at_j != tgtFmt.end(); ++at_j){
                float cut_j = at_j->pse->bondcut;
                if (cut_j<=0) {
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
    }

    void setBondsCell(float cutfac) const
    {
        if(getNat() >= 1000){
            setBondsCellSplit(cutfac);
        }else{
            setBondsCellTrivial(cutfac);
        }
    }
    void setBondsCellSplit(float cutfac) const
    {
        const auto asCrystal = asFmt(AtomFmt::Crystal);
        asCrystal.evaluateCache();
        enum wrapDir{x, y, z, xy, xmy, xz, xmz, yz, ymz, xyz, xymz, xmyz, mxyz};
        // get largest cutoff
        float cut{0};
        for (const auto& at:asCrystal) {
            cut = std::max(at.pse->bondcut, cut);
        }
        cut = 5*cut;
        auto cell = this->cell->cellvec * this->cell->dimBohr;
        // get cell lengths
        Vec size_split{};
        size_split[0] = std::abs(cell[0][0])+std::abs(cell[1][0])+std::abs(cell[2][0]);
        size_split[1] = std::abs(cell[0][1])+std::abs(cell[1][1])+std::abs(cell[2][1]);
        size_split[2] = std::abs(cell[0][2])+std::abs(cell[1][2])+std::abs(cell[2][2]);
        // fragment cell
        SizeVec n_split{1,1,1};
        if(size_split[0] >= cut){
            n_split[0] = std::max(size_t{1}, static_cast<size_t>(std::round(size_split[0] / cut)));
            size_split[0] = 1.f/n_split[0];
        }
        if(size_split[1] >= cut){
            n_split[1] = std::max(size_t{1}, static_cast<size_t>(std::round(size_split[1] / cut)));
            size_split[1] = 1.f/n_split[1];
        }
        if(size_split[2] >= cut){
            n_split[2] = std::max(size_t{1}, static_cast<size_t>(std::round(size_split[2] / cut)));
            size_split[2] = 1.f/n_split[2];
        }
        // put atoms in boxes
        DataGrid3D<std::vector<size_t>> boxes{n_split};
        for(auto it=asCrystal.begin(); it!=asCrystal.end(); ++it){
            auto findBox = [&](size_t dir){
                auto tmp = std::fmod(it->coord[dir], 1);
                if (tmp<0) tmp+=1;
                return static_cast<size_t>(tmp/size_split[dir]);
            };
            boxes(findBox(0), findBox(1), findBox(2)).push_back(it.getIdx());
        }
        // get bonds by iterating over boxes and their neighbors
        auto& bonds = this->bonds->bonds;
        // indices of neighbor boxes
        size_t xl, xh, yl, yh, zl, zh;
        // wrapping needed for specific neighbor
        bool w_xl{}, w_xh{}, w_yl{}, w_yh{}, w_zl{}, w_zh;
        for(size_t x=0; x<boxes.extent[0]; ++x){
            if(x == boxes.extent[0]-1){
                xh = 0; w_xh = true;
            }else{
                xh = x+1; w_xh = false;
            }
            if(x == 0){
                xl = boxes.extent[0]-1; w_xl = true;
            }else{
                xl = x-1; w_xl = false;
            }
            for(size_t y=0; y<boxes.extent[1]; ++y){
                if(y == boxes.extent[1]-1){
                    yh = 0; w_yh = true;
                }else{
                    yh = y+1; w_yh = false;
                }
                if(y == 0){
                    yl = boxes.extent[1]-1; w_yl = true;
                }else{
                    yl = y-1; w_yl = false;
                }
                for(size_t z=0; z<boxes.extent[2]; ++z){
                    if(z == boxes.extent[2]-1){
                        zh = 0; w_zh = true;
                    }else{
                        zh = z+1; w_zh = false;
                    }
                    if(z == 0){
                        zl = boxes.extent[2]-1; w_zl = true;
                    }else{
                        zl = z-1; w_zl = false;
                    }
                    for(auto it_i=boxes(x,y,z).begin(); it_i != boxes(x,y,z).end(); ++it_i){
                        auto i = *it_i;
                        const auto& at_i = asCrystal[i];
                        auto cut_i = at_i.pse->bondcut;
                        if (cut_i <= 0){
                            continue;
                        }
                        // loop over rest of current box
                        for (auto it_j = it_i+1; it_j!=boxes(x,y,z).end(); ++it_j) {
                            auto j = *it_j;
                            const auto& at_j = asCrystal[j];
                            auto cut_j = at_j.pse->bondcut;
                            if (cut_j <= 0) {
                                continue;
                            }
                            auto effcut = (cut_i + cut_j) * cutfac;
                            Vec dist_v = at_i.coord - at_j.coord;
                            std::array<int16_t, 3> diff_v, crit_v;
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
                            if(!(crit_v[0]|crit_v[1]|crit_v[2])){
                                // TODO: fail here? set flag? overlapping atoms!
                                continue;
                            }
                            // convert dist_v to bohr, wrap into cell if needed
                            for(size_t d=0; d<3; ++d){
                                if(std::abs(dist_v[d]) > 0.5f){
                                    dist_v[d] -= crit_v[d];
                                    diff_v[d] += crit_v[d];
                                }
                            }
                            dist_v = dist_v * cell;
                            if ((dist_v[0] > effcut) ||
                                (dist_v[1] > effcut) ||
                                (dist_v[2] > effcut)) {
                                continue;
                            }
                            auto dist_n = Vec_dot(dist_v, dist_v);
                            if((0.57f < dist_n) && (dist_n < effcut*effcut)) {
                                bonds.push_back({i, j, std::sqrt(dist_n),
                                                 diff_v[0], diff_v[1], diff_v[2]});
                            }
                        }
                        // other boxes should be evaluated in total
                        auto wrapcheck = [&](size_t xt, size_t yt, size_t zt){
                            for (const auto& j: boxes(xt, yt, zt)) {
                                const auto& at_j = asCrystal[j];
                                auto cut_j = at_j.pse->bondcut;
                                if (cut_j <= 0) {
                                    continue;
                                }
                                auto effcut = (cut_i + cut_j) * cutfac;
                                Vec dist_v = at_i.coord - at_j.coord;
                                std::array<int16_t, 3> diff_v, crit_v;
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
                                // convert dist_v to bohr, wrap into cell if needed
                                for(size_t d=0; d<3; ++d){
                                    if(std::abs(dist_v[d]) > 0.5f){
                                        dist_v[d] -= crit_v[d];
                                        diff_v[d] += crit_v[d];
                                    }
                                }
                                dist_v = dist_v * cell;
                                if ((dist_v[0] > effcut) ||
                                    (dist_v[1] > effcut) ||
                                    (dist_v[2] > effcut)) {
                                    continue;
                                }
                                auto dist_n = Vec_dot(dist_v, dist_v);
                                if((0.57f < dist_n) && (dist_n < effcut*effcut)) {
                                    bonds.push_back({i, j, std::sqrt(dist_n),
                                                     diff_v[0], diff_v[1], diff_v[2]});
                                }
                            }
                        };
                        // x
                        wrapcheck(xh, y, z);
                        // xy
                        wrapcheck(xh, yh, z);
                        // xyz
                        wrapcheck(xh, yh, zh);
                        // xy-z
                        wrapcheck(xh, yh, zl);
                        // x-y
                        wrapcheck(xh, yl, z);
                        // x-yz
                        wrapcheck(xh, yl, zh);
                        // xz
                        wrapcheck(xh, y, zh);
                        // x-z
                        wrapcheck(xh, y, zl);
                        // y
                        wrapcheck(x, yh, z);
                        // yz
                        wrapcheck(x, yh, zh);
                        // -xyz
                        wrapcheck(xl, yh, zh);
                        // y-z
                        wrapcheck(x, yh, zl);
                        // z
                        wrapcheck(x, y, zh);
                    }
                }
            }
        }
    }
    void setBondsCellTrivial(float cutfac) const
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
                // evaluation-lambda
                auto checkBond = [&](const Vec& dist, const std::array<int16_t, 3>& offset)
                {
                    if ((dist[0]>effcut) || (dist[1]>effcut) || (dist[2]>effcut)) {
                        return;
                    }
                    float dist_n = Vec_dot(dist, dist);
                    if ((0.57f < dist_n) && (dist_n < effcut*effcut)) {
                        bonds.push_back({i, j, std::sqrt(dist_n), offset[0], offset[1], offset[2]});
                    }
                };
                // 0-vector
                checkBond(dist_v, diff_v);
                if(crit_v[0] != 0){
                    // x, -x
                    checkBond(dist_v-crit_v[0]*x,
                              {{static_cast<int16_t>(diff_v[0]+crit_v[0]),diff_v[1],diff_v[2]}});
                }
                if(crit_v[1] != 0){
                    // y, -y
                    checkBond(dist_v-crit_v[1]*y,
                              {{diff_v[0],static_cast<int16_t>(diff_v[1]+crit_v[1]),diff_v[2]}});
                    if(crit_v[0] == crit_v[1]){
                        // x+y, -x-y
                        checkBond(dist_v-crit_v[0]*xy,
                                  {{static_cast<int16_t>(diff_v[0]+crit_v[0]),
                                    static_cast<int16_t>(diff_v[1]+crit_v[1]),
                                    diff_v[2]}});
                    }else if(crit_v[0] == -crit_v[1]){
                        // x-y, -x+y
                        checkBond(dist_v-crit_v[0]*xmy,
                                  {{static_cast<int16_t>(diff_v[0]+crit_v[0]),
                                    static_cast<int16_t>(diff_v[1]+crit_v[1]),
                                    diff_v[2]}});
                    }
                }
                if(crit_v[2] != 0){
                    // z, -z
                    checkBond(dist_v-crit_v[2]*z,
                              {{diff_v[0],diff_v[1],static_cast<int16_t>(diff_v[2]+crit_v[2])}});
                    if(crit_v[0] == crit_v[2]){
                        // x+z, -x-z
                        checkBond(dist_v-crit_v[0]*xz,
                                  {{static_cast<int16_t>(diff_v[0]+crit_v[0]),
                                    diff_v[1],
                                    static_cast<int16_t>(diff_v[2]+crit_v[2])}});
                    }else if(crit_v[0] == -crit_v[2]){
                        // x-z, -x+z
                        checkBond(dist_v-crit_v[0]*xmz,
                                  {{static_cast<int16_t>(diff_v[0]+crit_v[0]),
                                    diff_v[1],
                                    static_cast<int16_t>(diff_v[2]+crit_v[2])}});
                    }
                    if(crit_v[1] == crit_v[2]){
                        // y+z, -y-z
                        checkBond(dist_v-crit_v[1]*yz,
                                  {{diff_v[0],
                                    static_cast<int16_t>(diff_v[1]+crit_v[1]),
                                    static_cast<int16_t>(diff_v[2]+crit_v[2])}});
                        if(crit_v[0] == crit_v[2]){
                            // x+y+z, -x-y-z
                            checkBond(dist_v-crit_v[0]*xyz,
                                      {{static_cast<int16_t>(diff_v[0]+crit_v[0]),
                                        static_cast<int16_t>(diff_v[1]+crit_v[1]),
                                        static_cast<int16_t>(diff_v[2]+crit_v[2])}});
                        } else if(crit_v[0] == -crit_v[2]){
                            // x-y-z, -x+y+z
                            checkBond(dist_v-crit_v[1]*mxyz,
                                      {{static_cast<int16_t>(diff_v[0]+crit_v[0]),
                                        static_cast<int16_t>(diff_v[1]+crit_v[1]),
                                        static_cast<int16_t>(diff_v[2]+crit_v[2])}});
                        }
                    }else if(crit_v[1] == -crit_v[2]){
                        // y-z, -y+z
                        checkBond(dist_v-crit_v[1]*ymz,
                                  {{diff_v[0],
                                    static_cast<int16_t>(diff_v[1]+crit_v[1]),
                                    static_cast<int16_t>(diff_v[2]+crit_v[2])}});
                        if(crit_v[0] == crit_v[2]){
                            // x-y+z, -x+y-z
                            checkBond(dist_v-crit_v[0]*xmyz,
                                      {{static_cast<int16_t>(diff_v[0]+crit_v[0]),
                                        static_cast<int16_t>(diff_v[1]+crit_v[1]),
                                        static_cast<int16_t>(diff_v[2]+crit_v[2])}});
                        } else if(crit_v[0] == -crit_v[2]){
                            // x+y-z, -x-y+z
                            checkBond(dist_v-crit_v[0]*xymz,
                                      {{static_cast<int16_t>(diff_v[0]+crit_v[0]),
                                        static_cast<int16_t>(diff_v[1]+crit_v[1]),
                                        static_cast<int16_t>(diff_v[2]+crit_v[2])}});
                        }
                    }
                }
            }
        }
    }
};

}

#endif // STEPCONST_H
