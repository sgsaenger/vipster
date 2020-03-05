#ifndef STEPBASE_H
#define STEPBASE_H

#include "atom.h"
#include "bond.h"
#include "vec.h"
#include "data.h"
#include "stepformatter.h"
#include "stepselection.h"

#include <memory>
#include <vector>
#include <set>

namespace Vipster {

namespace detail{
/*
 * Helper functions to ensure step Formatters or Selections aren't nested in itself, e.g.:
 *
 * make_formatter_t<Formatter<T>> will return Formatter<T> instead of Formatter<Formatter<T>>
 *
 * elso enforce that Selection has priority over Formatter, i.e. Selection<Formatter<T>>
 */

template<typename, template<typename...> typename>
static constexpr bool is_instance{false};

template<typename A, template<typename...> typename B>
static constexpr bool is_instance<B<A>, B>{true};

template<typename T>
static constexpr bool is_selection = is_instance<T, Selection>;

template<typename T>
static constexpr bool is_formatter = is_instance<T, Formatter>;

template<typename T>
struct make_selection {
    using type = Selection<T>;
};

template<typename T>
struct make_selection<Selection<T>> {
    using type = typename make_selection<T>::type;
};

template<typename T>
using make_selection_t = typename make_selection<T>::type;

template<typename T>
struct make_formatter {
    using type = Formatter<T>;
};

template<typename T>
struct make_formatter<Formatter<T>>{
    using type = typename make_formatter<T>::type;
};

template<typename T>
struct make_formatter<Selection<T>>{
    using type = make_selection_t<typename make_formatter<T>::type>;
};

template<typename T>
using make_formatter_t = typename make_formatter<T>::type;

}

/*
 * Const-Base for all Step-like containers
 *
 * Implements const-interface for common state:
 * - Atoms (atomcontainer must be provided as template argument)
 * - Bonds
 * - Comment
 */
template<typename T>
class StepConst
{
    template<typename U> friend class StepConst;
public:
    virtual ~StepConst() = default;
    StepConst(const StepConst&) = default;
    StepConst(StepConst &&) = default;
    StepConst& operator=(StepConst&&) = default;
    StepConst& operator=(const StepConst&) = default;

    // Periodic table
    const PeriodicTable& getPTE() const
    {
        return *atoms->ctxt.pte;
    }

    // Atoms
    using atom_source = T;
    using const_atom = typename T::const_atom;
    size_t              getNat() const noexcept
    {
        return atoms->getNat();
    }
    const_atom          operator[](size_t i) const
    {
        return {*atoms, i};
    }
    const_atom          at(size_t i) const
    {
        if(i>=getNat()){
            throw Error("Atom-index out of bounds");
        }else{
            return operator[](i);
        }
    }
    const_atom          front() const
    {
        return {*atoms, 0};
    }
    const_atom          back() const
    {
        return {*atoms, getNat()-1};
    }
    const atom_source&  getAtoms() const noexcept
    {
        return *atoms;
    }
    // Atom-iterators
    using const_iterator = typename T::const_iterator;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;
    const_iterator          begin() const
    {
        return const_iterator{*atoms, 0};
    }
    const_iterator          cbegin() const
    {
        return const_iterator{*atoms, 0};
    }
    const_iterator          end() const
    {
        return const_iterator{*atoms, getNat()};
    }
    const_iterator          cend() const
    {
        return const_iterator{*atoms, getNat()};
    }
    const_reverse_iterator  rbegin() const
    {
        return std::make_reverse_iterator(end());
    }
    const_reverse_iterator  crbegin() const
    {
        return std::make_reverse_iterator(cend());
    }
    const_reverse_iterator  rend() const
    {
        return std::make_reverse_iterator(begin());
    }
    const_reverse_iterator  crend() const
    {
        return std::make_reverse_iterator(cbegin());
    }

    // Selection
    using const_selection = StepConst<detail::make_selection<T>>;
    const_selection select(SelectionFilter filter) const
    {
        if constexpr(detail::is_selection<T>){
            return {std::make_shared<typename const_selection::atom_source>(
                            this->atoms->atoms, this->atoms->indices),
                    this->bonds, this->comment};
        }else{
            return {std::make_shared<typename const_selection::atom_source>(
                            this->atoms, evalFilter(*this, filter)),
                    this->bonds, this->comment};
        }
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
    using const_formatter = StepConst<detail::make_formatter_t<T>>;
    const_formatter     asFmt(AtomFmt tgt) const
    {
        // create Formatter<T> from Formatter<T>
        if constexpr(detail::is_formatter<T>){
            return {std::make_shared<typename const_formatter::atom_source>(
                            this->atoms->atoms, tgt),
                    this->bonds, this->comment};
        }else if constexpr(detail::is_selection<T>){
            using base_source = std::remove_reference_t<decltype(*this->atoms->atoms)>;
            if constexpr(detail::is_formatter<base_source>){
                // create Selection<Formatter<T>> from Formatter<T>
                return {std::make_shared<typename const_formatter::atom_source>(
                            std::make_shared<base_source>(
                                this->atoms->atoms->atoms, tgt),
                            this->atoms->indices),
                        this->bonds, this->comment};
            }else{ // create Selection<Formatter<T>> from T
                return {std::make_shared<typename const_formatter::atom_source>(
                            std::make_shared<detail::Formatter<base_source>>(
                                this->atoms->atoms, tgt),
                            this->atoms->indices),
                        this->bonds, this->comment};
            }
        }else{ // create Formatter<T> from T
            return {std::make_shared<typename const_formatter::atom_source>(
                            this->atoms, tgt),
                    this->bonds, this->comment};
        }
    }
    AtomFmt             getFmt() const
    {
        return atoms->ctxt.fmt;
    }
    void                setFmt(AtomFmt tgt) const
    {
        if(atomFmtRelative(tgt)){
            atoms->ctxt.cell->enabled = true;
        }
        atoms->ctxt.fmt = tgt;
    }

    // Bonds
    const std::vector<Bond>&    getBonds() const
    {
        return bonds->bonds;
    }

    // Cell
    bool    hasCell() const
    {
        return atoms->ctxt.cell->enabled;
    }
    double  getCellDim(AtomFmt fmt) const
    {
        if(atomFmtRelative(fmt) || fmt > detail::AtomContext::toAngstrom.size()){
            throw Error{"StepConst::getCellDim: Invalid AtomFmt, needs to be absolute"};
        }
        return atoms->ctxt.cell->dimension * detail::AtomContext::fromAngstrom[fmt];
    }
    Mat     getCellVec() const
    {
        return atoms->ctxt.cell->matrix;
    }
    Vec     getCom() const
    {
        return getCom(atoms->ctxt.fmt);
    }
    Vec     getCom(AtomFmt fmt) const
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
        for(const auto& at: asFmt(fmt)){
            min[0]=std::min(min[0],at.coord[0]);
            min[1]=std::min(min[1],at.coord[1]);
            min[2]=std::min(min[2],at.coord[2]);
            max[0]=std::max(max[0],at.coord[0]);
            max[1]=std::max(max[1],at.coord[1]);
            max[2]=std::max(max[2],at.coord[2]);
        }
        return (min+max)/2;
    }
    Vec     getCenter(AtomFmt fmt, bool com=false) const
    {
        if(com || !atoms->ctxt.cell->enabled){
            return getCom(fmt);
        }
        const Mat& cv = atoms->ctxt.cell->matrix;
        return (cv[0]+cv[1]+cv[2]) * getCellDim(fmt) / 2;
    }

protected:
    StepConst(std::shared_ptr<atom_source> atoms, std::shared_ptr<BondList> bonds,
              std::shared_ptr<std::string> comment)
        : atoms{atoms}, bonds{bonds}, comment{comment}
    {}
    // Data
    std::shared_ptr<atom_source>    atoms;
    std::shared_ptr<BondList>       bonds;
    std::shared_ptr<std::string>    comment;
};

/*
 * Base for mutable Step-like containers
 *
 * Implements non-const interface and common utility functions
 *
 * Template-argument shall be atomcontainer
 */
template<typename T>
class StepMutable: public StepConst<T>
{
    template<typename U> friend class StepMutable;
public:
    virtual ~StepMutable() = default;
    StepMutable(const StepMutable &) = default;
    StepMutable(StepMutable &&) = default;
    StepMutable& operator=(const StepMutable&) = default;
    StepMutable& operator=(StepMutable&&) = default;

    // Periodic table
    PeriodicTable& getPTE()
    {
        return *this->atoms->ctxt.pte;
    }

    // Atoms
    using typename StepConst<T>::atom_source;
    using atom = typename T::atom;
    using StepConst<T>::operator[];
    atom        operator[](size_t i)
    {
        return {*this->atoms, i};
    }
    using StepConst<T>::at;
    atom        at(size_t i)
    {
        if(i >= this->getNat()){
            throw Error("Atom-index out of bounds");
        }else{
            return operator[](i);
        }
    }
    atom        front()
    {
        return {*this->atoms, 0};
    }
    atom        back()
    {
        return {*this->atoms, this->getNat()-1};
    }
    using StepConst<T>::begin;
    // Atom-iterators
    using iterator = typename T::iterator;
    using reverse_iterator = std::reverse_iterator<iterator>;
    iterator    begin()
    {
        return iterator{*this->atoms, 0};
    }
    using StepConst<T>::end;
    iterator    end()
    {
        return iterator{*this->atoms, this->getNat()};
    }
    using StepConst<T>::rbegin;
    reverse_iterator rbegin()
    {
        return std::make_reverse_iterator(end());
    }
    using StepConst<T>::rend;
    reverse_iterator rend()
    {
        return std::make_reverse_iterator(begin());
    }

    // Selection
    using StepConst<T>::select;
    using selection = StepMutable<detail::make_selection_t<T>>;
    selection select(SelectionFilter filter)
    {
        if constexpr(detail::is_selection<T>){
            return {std::make_shared<typename selection::atom_source>(this->atoms->atoms, this->atoms->indices),
                    this->bonds, this->comment};
        }else{
            return {std::make_shared<typename selection::atom_source>(this->atoms, evalFilter(*this, filter)),
                    this->bonds, this->comment};
        }
    }

    // Comment
    void                setComment(const std::string& s)
    {
        *this->comment = s;
    }

    // Format
    using StepConst<T>::asFmt;
    using formatter = StepMutable<detail::make_formatter_t<T>>;
    formatter     asFmt(AtomFmt tgt)
    {
        // create Formatter<T> from Formatter<T>
        if constexpr(detail::is_formatter<T>){
            return {std::make_shared<typename formatter::atom_source>(
                            this->atoms->atoms, tgt),
                    this->bonds, this->comment};
        }else if constexpr(detail::is_selection<T>){
            using base_source = std::remove_reference_t<decltype(*this->atoms->atoms)>;
            if constexpr(detail::is_formatter<base_source>){
                // create Selection<Formatter<T>> from Formatter<T>
                return {std::make_shared<typename formatter::atom_source>(
                            std::make_shared<base_source>(
                                this->atoms->atoms->atoms, tgt),
                            this->atoms->indices),
                        this->bonds, this->comment};
            }else{ // create Selection<Formatter<T>> from T
                return {std::make_shared<typename formatter::atom_source>(
                            std::make_shared<detail::Formatter<base_source>>(
                                this->atoms->atoms, tgt),
                            this->atoms->indices),
                        this->bonds, this->comment};
            }
        }else{ // create Formatter<T> from T
            return {std::make_shared<typename formatter::atom_source>(
                            this->atoms, tgt),
                    this->bonds, this->comment};
        }
    }

    // Bonds
    void setBonds() const
    {
        if(this->atoms->ctxt.cell->enabled){
            setBondsCell();
        }else{
            setBondsMolecule();
        }
    }
    void addBond(size_t at1, size_t at2, DiffVec diff={}, const std::string& type="")
    {
        // calc distance in bohr
        auto getDistance = [&](size_t at1, size_t at2, DiffVec diff)
        {
            auto fmt = asFmt(AtomFmt::Bohr);
            if(std::any_of(diff.begin(), diff.end(), [](auto i)->bool{return i;})){
                // have to calculate periodic bond
                auto cv = this->getCellVec() * this->getCellDim(AtomFmt::Bohr);
                return Vec_length(fmt[at1].coord - fmt[at2].coord
                                  - diff[0] * cv[0] - diff[1] * cv[1]
                                  - diff[2] * cv[2]);
            }else{
                // just use immediate distance
                return Vec_length(fmt[at1].coord - fmt[at2].coord);
            }
        };
        if(!type.empty()){
            // register/look-up type, then create bond
            this->bonds->bonds.push_back({at1, at2, getDistance(at1, at2, diff), diff,
                                            &*this->bonds->types.emplace(type,
                                             defaultColors[this->bonds->types.size()%5]).first});
        }else{
            // create untyped bond
            this->bonds->bonds.push_back({at1, at2, getDistance(at1, at2, diff), diff,
                                            nullptr});
        }
    }
    void delBond(size_t idx)
    {
        auto& bonds = this->bonds->bonds;
        bonds.erase(bonds.begin()+idx);
    }
    void setBondType(size_t idx, std::string type)
    {
        auto& bond = this->bonds->bonds[idx];
        if(type.empty()){
            bond.type = nullptr;
        }else{
            bond.type = &*this->bonds->types.emplace(
                        type,
                        defaultColors[this->bonds->types.size()%5]
                        ).first;
        }
    }

    // Modifier functions
    void modScale(AtomFmt tgt){
        if(tgt == this->atoms->ctxt.fmt){ return; }
        auto tmp = asFmt(tgt);
        auto source = tmp.cbegin();
        auto target = this->begin();
        while(source != tmp.cend()){
            target->coord = source->coord;
            ++target;
            ++source;
        }
        this->atoms->ctxt.fmt = tgt;
    }
    void modShift(Vec shift, double fac=1.0){
        shift *= fac;
        for(auto& at:*this){
            at.coord += shift;
        }
    }
    void modRotate(double angle, Vec axis, Vec shift={0,0,0}){
        angle *= deg2rad;
        double c = std::cos(angle);
        double s = -std::sin(angle);
        double ic = 1.-c;
        auto crystal = this->atoms->ctxt.fmt == AtomFmt::Crystal;
        auto copy = this->atoms->ctxt;
        copy.fmt = AtomFmt::Bohr;
        auto toBohr = detail::makeConverter(this->atoms->ctxt, copy);
        if(crystal) axis = toBohr(axis);
        double len = Vec_length(axis);
        if(float_comp(len, 0.)){
            throw Error("0-Vector cannot be rotation axis");
        }
        axis /= len;
        Mat rotMat = {Vec{ic * axis[0] * axis[0] + c,
                          ic * axis[0] * axis[1] - s * axis[2],
                          ic * axis[0] * axis[2] + s * axis[1]},
                      Vec{ic * axis[1] * axis[0] + s * axis[2],
                          ic * axis[1] * axis[1] + c,
                          ic * axis[1] * axis[2] - s * axis[0]},
                      Vec{ic * axis[2] * axis[0] - s * axis[1],
                          ic * axis[2] * axis[1] + s * axis[0],
                          ic * axis[2] * axis[2] + c}};
        if(!crystal){
            for(auto& at:*this){
                at.coord = (at.coord - shift) * rotMat + shift;
            }
        }else{
            const auto fromBohr = detail::makeConverter(copy, this->atoms->ctxt);
            shift = toBohr(shift);
            for(auto &at: *this){
                at.coord = fromBohr((toBohr(at.coord) - shift) * rotMat + shift);
            }
        }
    }
    void modMirror(Vec ax1, Vec ax2, Vec shift={0,0,0}){
        auto crystal = this->atoms->ctxt.fmt == AtomFmt::Crystal;
        auto copy = this->atoms->ctxt;
        copy.fmt = AtomFmt::Bohr;
        auto toBohr = detail::makeConverter(this->atoms->ctxt, copy);
        if(crystal){
            ax1 = toBohr(ax1);
            ax2 = toBohr(ax2);
            shift = toBohr(shift);
        }
        Vec normal = Vec_cross(ax1, ax2);
        normal /= Vec_length(normal);
        if(!crystal){
            for(auto& at:*this){
                at.coord -= 2*Vec_dot(at.coord-shift, normal)*normal;
            }
        }else{
            auto fromBohr = detail::makeConverter(copy, this->atoms->ctxt);
            for(auto &at: *this){
                at.coord -= fromBohr(2*Vec_dot(toBohr(at.coord) - shift, normal) * normal);
            }
        }
    }

protected:
    StepMutable(std::shared_ptr<T> atoms, std::shared_ptr<BondList> bonds,
                std::shared_ptr<std::string> comment)
        : StepConst<T>{atoms, bonds, comment}
    {}

private:

    // generate non-periodic bonds
    void setBondsMolecule() const
    {
        // get suitable absolute representation
        const AtomFmt fmt = (this->atoms->ctxt.fmt == AtomFmt::Angstrom) ? AtomFmt::Angstrom : AtomFmt::Bohr;
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
        if(float_comp(cut, 0.)){
            // if we have no bonding atom types,
            // exit early as no bonds will be found
            return;
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
        DataGrid<3, std::vector<size_t>> bins{n_split};
        if(n_split == SizeVec{1,1,1}){
            // only one bin
            auto& bin = bins(0,0,0);
            bin.resize(this->getNat());
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

    // generate periodic bonds
    void setBondsCell() const
    {
        const auto asCrystal = asFmt(AtomFmt::Crystal);
        auto& bonds = this->bonds->bonds;
        const auto cell = this->getCellVec() * this->getCellDim(AtomFmt::Bohr);
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
        DataGrid<3, std::vector<size_t>> bins{n_split};
        if(n_split == SizeVec{1,1,1}){
            // only one bin
            auto& bin = bins(0,0,0);
            bin.resize(this->getNat());
            std::iota(bin.begin(), bin.end(), 0);
        }else{
            // multiple bins
            for(auto it = asCrystal.begin(); it != asCrystal.end(); ++it){
                auto findBin = [&](size_t dir){
                    auto tmp = std::fmod(it->coord[dir], 1);
                    if (tmp<0) tmp+=1;
                    auto tgt = static_cast<size_t>(tmp/size_split[dir]);
                    if (tgt==n_split[dir]) tgt-=1;
                    return tgt;
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

#endif // STEPBASE_H
