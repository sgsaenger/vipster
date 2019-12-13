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
#include <optional>

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
    std::shared_ptr<PeriodicTable> pte;

    // Selection
    constSelection select(std::string filter) const
    {
        return {pte, at_fmt, this, filter, cell, comment};
    }
    constSelection select(SelectionFilter filter) const
    {
        return {pte, at_fmt, this, filter, cell, comment};
    }

    // Atoms
    size_t          getNat() const noexcept
    {
        return atoms->getNat();
    }
    constAtom       operator[](size_t i) const noexcept
    {
        return *const_iterator{atoms, pte, at_fmt, i};
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
        return const_iterator{atoms, pte, at_fmt, 0};
    }
    const_iterator   cbegin() const noexcept
    {
        return const_iterator{atoms, pte, at_fmt, 0};
    }
    const_iterator   end() const noexcept
    {
        return const_iterator{atoms, pte, at_fmt, getNat()};
    }
    const_iterator   cend() const noexcept
    {
        return const_iterator{atoms, pte, at_fmt, getNat()};
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
        // evaluate tgt-cache via side-effect
        asFmt(tgt);
        at_fmt = tgt;
    }
    StepConst           asFmt(AtomFmt tgt) const
    {
        auto tmp = StepConst{pte, tgt, atoms, bonds, cell, comment};
        tmp.evaluateCache();
        return tmp;
    }
    std::function<Vec(const Vec&)>  getFormatter(AtomFmt source,
                                                 AtomFmt target) const
    {
        double fac{};
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
    const std::vector<Bond>&    getBonds() const
    {
        // TODO: was passiert wenn bindungen manuell sind, aber atome gelöscht werden?
        if((bonds->mode == BondMode::Automatic) &&
            bonds->outdated)
        {
            setBonds();
        }
        return bonds->bonds;
    }
    void                        setBonds() const
    {
        bonds->bonds.clear();
        if(!getNat()){
            bonds->outdated = false;
            return;
        }
        if(cell->enabled){
            setBondsCell();
        }else{
            setBondsMolecule();
        }
        bonds->outdated = false;
    }
    void                        setBondMode(BondMode mode) const
    {
        bonds->mode = mode;
    }
    BondMode                    getBondMode() const
    {
        return bonds->mode;
    }

    // Cell
    bool    hasCell() const noexcept
    {
        return cell->enabled;
    }
    double   getCellDim(CdmFmt fmt) const noexcept
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
        Vec min{{std::numeric_limits<double>::max(),
                 std::numeric_limits<double>::max(),
                 std::numeric_limits<double>::max()}};
        Vec max{{std::numeric_limits<double>::lowest(),
                 std::numeric_limits<double>::lowest(),
                 std::numeric_limits<double>::lowest()}};
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
    StepConst(std::shared_ptr<PeriodicTable> pte, AtomFmt fmt,
             std::shared_ptr<source> atoms, std::shared_ptr<BondList> bonds,
             std::shared_ptr<CellData> cell, std::shared_ptr<std::string> comment)
        : pte{pte}, at_fmt{fmt}, atoms{atoms}, bonds{bonds}, cell{cell}, comment{comment}
    {}
    // Data
    mutable AtomFmt                 at_fmt; // mutable, only controls visibility
    std::shared_ptr<source>         atoms;  // immutable
    std::shared_ptr<BondList>       bonds;  // mutable, is only cache
    std::shared_ptr<CellData>       cell;   // immutable except enabled/disable-state
    std::shared_ptr<std::string>    comment; // immutable

private:

    // Bonds
    void setBondsMolecule() const
    {
        // get suitable absolute representation
        const AtomFmt fmt = (this->at_fmt == AtomFmt::Angstrom) ? AtomFmt::Angstrom : AtomFmt::Bohr;
        const double fmtscale{(fmt == AtomFmt::Angstrom) ? invbohr : 1};
        auto tgtFmt = asFmt(fmt);
        // get bounds of system and largest cutoff
        Vec min{0,0,0};
        Vec max{0,0,0};
        double cut{0};
        for (const auto& at:tgtFmt) {
            min[0] = std::min(at.coord[0], min[0]);
            min[1] = std::min(at.coord[1], min[1]);
            min[2] = std::min(at.coord[2], min[2]);
            max[0] = std::max(at.coord[0], max[0]);
            max[1] = std::max(at.coord[1], max[1]);
            max[2] = std::max(at.coord[2], max[2]);
            cut = std::max(at.type->bondcut, cut);
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
        // assign atoms to bins
        DataGrid3D<std::vector<size_t>> bins{n_split};
        if(n_split == SizeVec{1,1,1}){
            // only one bin
            auto& bin = bins(0,0,0);
            bin.resize(getNat());
            std::iota(bin.begin(), bin.end(), 0);
        }else{
            // multiple bins
            for(auto it = tgtFmt.begin(); it != tgtFmt.end(); ++it){
                auto findBin = [&](size_t dir){
                    return std::max(size_t{0},
                                    std::min(n_split[dir]-1,
                                             static_cast<size_t>((it->coord[dir]-min[dir])/size_split[dir])));
                };
                bins(findBin(0), findBin(1), findBin(2)).push_back(it.getIdx());
            }
        }
        // get bonds by iterating over bins and their neighbors
        auto& bonds = this->bonds->bonds;
        for(size_t x=0; x<bins.extent[0]; ++x){
            for(size_t y=0; y<bins.extent[1]; ++y){
                for(size_t z=0; z<bins.extent[2]; ++z){
                    for (auto it_i=bins(x,y,z).begin(); it_i != bins(x,y,z).end(); ++it_i) {
                        auto i = *it_i;
                        const auto& at_i = tgtFmt[i];
                        auto cut_i = at_i.type->bondcut;
                        if (cut_i <= 0){
                            continue;
                        }
                        // loop over rest of current bin
                        for (auto it_j = it_i+1; it_j!=bins(x,y,z).end(); ++it_j) {
                            auto j = *it_j;
                            const auto& at_j = tgtFmt[j];
                            auto cut_j = at_j.type->bondcut;
                            if (cut_j <= 0) {
                                continue;
                            }
                            auto effcut = (cut_i + cut_j) * 1.1f;
                            Vec dist_v = at_i.coord - at_j.coord;
                            if (((dist_v[0] *= fmtscale) > effcut) ||
                                ((dist_v[1] *= fmtscale) > effcut) ||
                                ((dist_v[2] *= fmtscale) > effcut)) {
                                continue;
                            }
                            auto dist_n = Vec_dot(dist_v, dist_v);
                            if((0.57f < dist_n) && (dist_n < effcut*effcut)) {
                                bonds.push_back({i, j, std::sqrt(dist_n), {}});
                            }
                        }
                        // visit neighboring bins
                        auto bondcheck = [&](size_t xt, size_t yt, size_t zt){
                            for (const auto& j: bins(xt, yt, zt)) {
                                const auto& at_j = tgtFmt[j];
                                auto cut_j = at_j.type->bondcut;
                                if (cut_j <= 0) {
                                    continue;
                                }
                                auto effcut = (cut_i + cut_j) * 1.1f;
                                Vec dist_v = at_i.coord - at_j.coord;
                                if (((dist_v[0] *= fmtscale) > effcut) ||
                                    ((dist_v[1] *= fmtscale) > effcut) ||
                                    ((dist_v[2] *= fmtscale) > effcut)) {
                                    continue;
                                }
                                auto dist_n = Vec_dot(dist_v, dist_v);
                                if((0.57f < dist_n) && (dist_n < effcut*effcut)) {
                                    bonds.push_back({i, j, std::sqrt(dist_n), {}});
                                }
                            }
                        };
                        // +x (includes -x)
                        if (x < bins.extent[0]-1) {
                            bondcheck(x+1, y, z);
                            // +x+y (includes -x-y)
                            if (y < bins.extent[1]-1) {
                                bondcheck(x+1, y+1, z);
                                // +x+y+z (includes -x-y-z)
                                if (z < bins.extent[2]-1) {
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
                                if (z < bins.extent[2]-1) {
                                    bondcheck(x+1, y-1, z+1);
                                }
                                // +x-y-z (includes -x+y+z)
                                if (z > 0) {
                                    bondcheck(x+1, y-1, z-1);
                                }
                            }
                            // +x+z (includes -x-z)
                            if (z < bins.extent[2]-1) {
                                bondcheck(x+1, y, z+1);
                            }
                            // +x-z (includex -x+z)
                            if (z > 0) {
                                bondcheck(x+1, y, z-1);
                            }
                        }
                        // +y (includes -y)
                        if (y < bins.extent[1]-1) {
                            bondcheck(x, y+1, z);
                            // +y+z (includes -y-z)
                            if (z < bins.extent[2]-1) {
                                bondcheck(x, y+1, z+1);
                            }
                            // +y-z (includes -y+z)
                            if (z > 0) {
                                bondcheck(x, y+1, z-1);
                            }
                        }
                        // +z (includes -z)
                        if (z < bins.extent[2]-1) {
                            bondcheck(x, y, z+1);
                        }
                    }
                }
            }
        }
    }

    void setBondsCell() const
    {
        const auto asCrystal = asFmt(AtomFmt::Crystal);
        auto& bonds = this->bonds->bonds;
        const auto cell = getCellVec() * getCellDim(CdmFmt::Bohr);
        // offset vectors for bin-bin interactions
        const Vec x_v = cell[0];
        const Vec y_v = cell[1];
        const Vec z_v = cell[2];
        const Vec xy_v   = x_v+y_v;
        const Vec xmy_v  = x_v-y_v;
        const Vec xz_v   = x_v+z_v;
        const Vec xmz_v  = x_v-z_v;
        const Vec yz_v   = y_v+z_v;
        const Vec ymz_v  = y_v-z_v;
        const Vec xyz_v  = xy_v + z_v;
        const Vec xymz_v = xy_v - z_v;
        const Vec xmyz_v = xz_v - y_v;
        const Vec mxyz_v = yz_v - x_v;
        // get largest cutoff times five, as criterion for partitioning the cell
        double cut{0};
        for(const auto& at: asCrystal) {
            cut = std::max(at.type->bondcut, cut);
        }
        cut = 5*cut;
        // get cell lengths
        Vec size_split{
            std::abs(cell[0][0]) + std::abs(cell[1][0]) + std::abs(cell[2][0]),
            std::abs(cell[0][1]) + std::abs(cell[1][1]) + std::abs(cell[2][1]),
            std::abs(cell[0][2]) + std::abs(cell[1][2]) + std::abs(cell[2][2])
        };
        // partition the cell
        SizeVec n_split{1, 1, 1};
        for(size_t i{0}; i<3; ++i){
            if(size_split[i] >= cut){
                n_split[i] = std::max(size_t{1}, static_cast<size_t>(std::round(size_split[i]/cut)));
                size_split[i] = 1./n_split[i];
            }
        }
        // assign atoms to bins
        DataGrid3D<std::vector<size_t>> bins{n_split};
        if(n_split == SizeVec{1,1,1}){
            // only one bin
            auto& bin = bins(0,0,0);
            bin.resize(getNat());
            std::iota(bin.begin(), bin.end(), 0);
        }else{
            // multiple bins
            for(auto it = asCrystal.begin(); it != asCrystal.end(); ++it){
                auto findBin = [&](size_t dir){
                    auto tmp = std::fmod(it->coord[dir], 1);
                    if (tmp<0) tmp+=1;
                    return static_cast<size_t>(tmp/size_split[dir]);
                };
                bins(findBin(0), findBin(1), findBin(2)).push_back(it.getIdx());
            }
        }
        /* we're going to loop over all neighboring bins,
         * unless any direction has only two bins
         */
        bool dir_periodic[3] {n_split[0] > 2, n_split[1] > 2, n_split[2] > 2};
        /* if any direction has only one bin,
         * we need to calculate (some) periodic self-interactions of this bin
         */
        bool dir_with_self[3] {n_split[0] == 1, n_split[1] == 1, n_split[2] == 1};
        /* if all directions have multiple bins,
         * periodic bonds exist only between bins,
         * which allows for some shortcuts
         */
        bool bin_with_self = dir_with_self[0] || dir_with_self[1] || dir_with_self[2];
        // iterate over source-bins
        for(size_t x=0; x<n_split[0]; ++x){
            size_t xh = (x < (n_split[0]-1)) ? x+1 : 0;
            size_t xl = (x > 0) ? x-1 : n_split[0]-1;
            for(size_t y=0; y<n_split[1]; ++y){
                size_t yh = (y < (n_split[1]-1)) ? y+1 : 0;
                size_t yl = (y > 0) ? y-1 : n_split[1]-1;
                for(size_t z=0; z<n_split[2]; ++z){
                    size_t zh = (z < (n_split[2]-1)) ? z+1 : 0;
                    size_t zl = (z > 0) ? z-1 : n_split[2]-1;
                    const auto& bin_source = bins(x,y,z);
                    for(auto it_source = bin_source.begin(); it_source != bin_source.end(); ++it_source){
                        auto idx_source = *it_source;
                        const auto& at_source = asCrystal[idx_source];
                        auto cut_source = at_source.type->bondcut;
                        if (cut_source <= 0){
                            // non-bonding atom type
                            continue;
                        }
                        // evaluate bonds between bin_source and bin(x_t,y_t,z_t)
                        auto checkBin = [&](size_t x_t, size_t y_t, size_t z_t){
                            const auto& bin_target = bins(x_t, y_t, z_t);
                            /* if bin_target == bin_source,
                             * we only need to visit atoms that have not been visited by the outer loop
                             * else we need to loop over all atoms in bin_target
                             */
                            auto it_target = (x_t == x && y_t == y && z_t == z) ?
                                        it_source+1 :
                                        bin_target.begin();
                            for(; it_target != bin_target.end(); ++it_target){
                                auto idx_target = *it_target;
                                const auto& at_target = asCrystal[idx_target];
                                auto cut_target = at_target.type->bondcut;
                                if (cut_target <= 0) {
                                    // non-bonding atom type
                                    continue;
                                }
                                // effective cutoff, with a 10% stretching margin
                                auto effcut = (cut_source + cut_target) * 1.1f;
                                // distance in crystal-units
                                Vec dist_v = at_source.coord - at_target.coord;
                                // diff_v contains whole-cell distance (truncated by casting to int)
                                DiffVec diff_v{
                                    static_cast<DiffVec::value_type>(dist_v[0]),
                                    static_cast<DiffVec::value_type>(dist_v[1]),
                                    static_cast<DiffVec::value_type>(dist_v[2])};
                                // wrap dist_v inside one cell
                                dist_v[0] = std::fmod(dist_v[0], 1);
                                dist_v[1] = std::fmod(dist_v[1], 1);
                                dist_v[2] = std::fmod(dist_v[2], 1);
                                // crit_v contains direction of distance-vector
                                auto getSignFuzzy = [](double f){
                                    return (std::abs(f) < std::numeric_limits<double>::epsilon()) ?
                                                0 : ((f<0) ? -1 : 1);
                                };
                                DiffVec crit_v{
                                    getSignFuzzy(dist_v[0]),
                                    getSignFuzzy(dist_v[1]),
                                    getSignFuzzy(dist_v[2]),
                                };
                                if(!(crit_v[0]|crit_v[1]|crit_v[2])){
                                    // TODO: fail here? set flag? overlapping atoms!
                                    continue;
                                }
                                // evaluation-lambda
                                auto checkBond = [&](const Vec& dist, const DiffVec& offset)
                                {
                                    if ((dist[0]>effcut) || (dist[1]>effcut) || (dist[2]>effcut)) {
                                        return;
                                    }
                                    double dist_n = Vec_dot(dist, dist);
                                    if ((0.57 < dist_n) && (dist_n < effcut*effcut)) {
                                        bonds.push_back({idx_source, idx_target, std::sqrt(dist_n), offset});
                                    }
                                };
                                /* if there are multiple bins in all directions
                                 * we can perform a naive wrapping,
                                 * as no bin-bin bonds with ourself are possible
                                 */
                                if (!bin_with_self) {
                                    for (size_t i=0; i<3; ++i){
                                        if(std::abs(dist_v[i]) > 0.5){
                                            dist_v[i] -= crit_v[i];
                                            diff_v[i] += crit_v[i];
                                        }
                                    }
                                }
                                // convert to bohr
                                dist_v = dist_v * cell;
                                // inside bin
                                checkBond(dist_v, diff_v);
                                // shortcut if no bin-bin interactions are possible
                                if (!bin_with_self){
                                    continue;
                                }
                                // x, -x
                                if (dir_with_self[0]) {
                                    if (crit_v[0] != 0) {
                                        checkBond(dist_v-crit_v[0]*x_v,
                                                  {{diff_v[0]+crit_v[0],diff_v[1],diff_v[2]}});
                                    }
                                }
                                // y, -y
                                if (dir_with_self[1]) {
                                    if (crit_v[1] != 0) {
                                        checkBond(dist_v-crit_v[1]*y_v,
                                                  {{diff_v[0],diff_v[1]+crit_v[1],diff_v[2]}});
                                        if (dir_with_self[0]) {
                                            // x+y, -x-y
                                            if (crit_v[0] == crit_v[1]) {
                                                checkBond(dist_v-crit_v[0]*xy_v,
                                                          {{diff_v[0]+crit_v[0],
                                                            diff_v[1]+crit_v[1],
                                                            diff_v[2]}});
                                            // x-y, -x+y
                                            } else if (crit_v[0] == -crit_v[1]) {
                                                checkBond(dist_v-crit_v[0]*xmy_v,
                                                          {{diff_v[0]+crit_v[0],
                                                            diff_v[1]+crit_v[1],
                                                            diff_v[2]}});
                                            }
                                        }
                                    }
                                }
                                // z, -z
                                if (dir_with_self[2]) {
                                    if (crit_v[2] != 0) {
                                        checkBond(dist_v-crit_v[2]*z_v,
                                                  {{diff_v[0],diff_v[1],diff_v[2]+crit_v[2]}});
                                        if (dir_with_self[0]) {
                                            // x+z, -x-z
                                            if (crit_v[0] == crit_v[2]) {
                                                checkBond(dist_v-crit_v[0]*xz_v,
                                                          {{diff_v[0]+crit_v[0],
                                                            diff_v[1],
                                                            diff_v[2]+crit_v[2]}});
                                            // x-z, -x+z
                                            } else if (crit_v[0] == -crit_v[2]) {
                                                checkBond(dist_v-crit_v[0]*xmz_v,
                                                          {{diff_v[0]+crit_v[0],
                                                            diff_v[1],
                                                            diff_v[2]+crit_v[2]}});
                                            }
                                        }
                                        if (dir_with_self[1]) {
                                            // y+z, -y-z
                                            if (crit_v[1] == crit_v[2]) {
                                                checkBond(dist_v-crit_v[1]*yz_v,
                                                          {{diff_v[0],
                                                            diff_v[1]+crit_v[1],
                                                            diff_v[2]+crit_v[2]}});
                                                if (dir_with_self[0]) {
                                                    // x+y+z, -x-y-z
                                                    if (crit_v[0] == crit_v[2]) {
                                                        checkBond(dist_v-crit_v[0]*xyz_v,
                                                                  {{diff_v[0]+crit_v[0],
                                                                    diff_v[1]+crit_v[1],
                                                                    diff_v[2]+crit_v[2]}});
                                                    // x-y-z, -x+y+z
                                                    } else if (crit_v[0] == -crit_v[2]) {
                                                        checkBond(dist_v-crit_v[1]*mxyz_v,
                                                                  {{diff_v[0]+crit_v[0],
                                                                    diff_v[1]+crit_v[1],
                                                                    diff_v[2]+crit_v[2]}});
                                                    }
                                                }
                                            // y-z, -y+z
                                            } else if (crit_v[1] == -crit_v[2]) {
                                                checkBond(dist_v-crit_v[1]*ymz_v,
                                                          {{diff_v[0],
                                                            diff_v[1]+crit_v[1],
                                                            diff_v[2]+crit_v[2]}});
                                                if (dir_with_self[0]) {
                                                    // x-y+z, -x+y-z
                                                    if (crit_v[0] == crit_v[2]) {
                                                        checkBond(dist_v-crit_v[0]*xmyz_v,
                                                                  {{diff_v[0]+crit_v[0],
                                                                    diff_v[1]+crit_v[1],
                                                                    diff_v[2]+crit_v[2]}});
                                                    // x+y-z, -x-y+z
                                                    } else if (crit_v[0] == -crit_v[2]) {
                                                        checkBond(dist_v-crit_v[0]*xymz_v,
                                                                  {{diff_v[0]+crit_v[0],
                                                                    diff_v[1]+crit_v[1],
                                                                    diff_v[2]+crit_v[2]}});
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        };
                        // visit current bin
                        checkBin(x,y,z);
                        // x
                        if (dir_periodic[0] || (x < (n_split[0]-1))) {
                            checkBin(xh, y, z);
                            // xy
                            if (dir_periodic[1] || (y < (n_split[1]-1))) {
                                checkBin(xh, yh, z);
                                // xyz
                                if (dir_periodic[2] || (z < (n_split[2]-1))) {
                                    checkBin(xh, yh, zh);
                                }
                                // xy-z
                                if (dir_periodic[2] || (z > 0)) {
                                    checkBin(xh, yh, zl);
                                }
                            }
                            // x-y
                            if (dir_periodic[1] || (y > 0)) {
                                checkBin(xh, yl, z);
                                // x-yz
                                if (dir_periodic[2] || (z < (n_split[2]-1))) {
                                    checkBin(xh, yl, zh);
                                }
                            }
                            // xz
                            if (dir_periodic[2] || (z < (n_split[2]-1))) {
                                checkBin(xh, y, zh);
                            }
                            // xy-z
                            if (dir_periodic[2] || (z > 0)) {
                                checkBin(xh, y, zl);
                            }
                        }
                        // y
                        if (dir_periodic[1] || (y < (n_split[1]-1))) {
                            checkBin(x, yh, z);
                            // yz
                            if (dir_periodic[2] || (z < (n_split[2]-1))) {
                                checkBin(x, yh, zh);
                                if (dir_periodic[0] || (x > 0)) {
                                    checkBin(xl, yh, zh);
                                }
                            }
                            // y-z
                            if (dir_periodic[2] || (z > 0)) {
                                checkBin(x, yh, zl);
                            }
                        }
                        // z
                        if (dir_periodic[2] || (z < (n_split[2]-1))) {
                            checkBin(x, y, zh);
                        }
                    }
                }
            }
        }
    }
};

}

#endif // STEPCONST_H
