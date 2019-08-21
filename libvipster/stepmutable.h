#ifndef STEPMUTABLE_H
#define STEPMUTABLE_H

#include "stepconst.h"

namespace Vipster {

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
public:
    virtual ~StepMutable() = default;

    using iterator = typename T::iterator;
    using reverse_iterator = std::reverse_iterator<iterator>;
    using selection = SelectionBase<StepMutable, T>;

    // Selection
    selection select(std::string filter)
    {
        return selection{this->pte, this->at_fmt, this, filter, this->cell, this->comment};
    }
    selection select(SelectionFilter filter)
    {
        return selection{this->pte, this->at_fmt, this, filter, this->cell, this->comment};
    }

    // Comment
    void                setComment(const std::string& s)
    {
        *this->comment = s;
    }

    // Format
    using StepConst<T>::asFmt;
    StepMutable asFmt(AtomFmt tgt)
    {
        auto tmp = StepMutable{this->pte, tgt,
                               this->atoms, this->bonds,
                               this->cell, this->comment};
        tmp.evaluateCache();
        return tmp;
    }

    // Atoms
    using StepConst<T>::operator[];
    Atom        operator[](size_t i) noexcept
    {
        return *iterator{this->atoms, this->pte, this->at_fmt, i};
    }
    using StepConst<T>::at;
    Atom        at(size_t i)
    {
        if(i >= this->getNat()){
            throw Error("Atom-index out of bounds");
        }else{
            return operator[](i);
        }
    }
    using StepConst<T>::begin;
    iterator    begin() noexcept
    {
        return iterator{this->atoms, this->pte, this->at_fmt, 0};
    }
    using StepConst<T>::end;
    iterator    end() noexcept
    {
        return iterator{this->atoms, this->pte, this->at_fmt, this->getNat()};
    }
    using StepConst<T>::rbegin;
    reverse_iterator rbegin() noexcept
    {
        return std::make_reverse_iterator(end());
    }
    using StepConst<T>::rend;
    reverse_iterator rend() noexcept
    {
        return std::make_reverse_iterator(begin());
    }

    // Bonds
    void newBond(size_t at1, size_t at2, DiffVec diff={}, const std::string* type=nullptr)
    {
        // calc distance in bohr
        auto getDistance = [&](size_t at1, size_t at2, DiffVec diff)
        {
            if(std::any_of(diff.begin(), diff.end(), [](auto i)->bool{return i;})){
                // have to calculate periodic bond
                auto cv = [&](){
                    switch(this->at_fmt){
                    case AtomFmt::Crystal:
                        return Mat{{{{1,0,0}}, {{0,1,0}}, {{0,0,1}}}};
                    case AtomFmt::Alat:
                        return this->getCellVec();
                    case AtomFmt::Angstrom:
                        return this->getCellVec() * this->getCellDim(CdmFmt::Angstrom);
                    case AtomFmt::Bohr:
                        return this->getCellVec() * this->getCellDim(CdmFmt::Bohr);
                    }
                    throw Error("Invalid AtomFmt");
                }();
                return Vec_length(this->formatVec(
                        operator[](at1).coord - operator[](at2).coord
                        - diff[0]*cv[0] - diff[1]*cv[1] - diff[2]*cv[2],
                        this->at_fmt, AtomFmt::Bohr));
            }else{
                // just use immediate distance
                return Vec_length(this->formatVec(
                        operator[](at1).coord - operator[](at2).coord,
                        this->at_fmt, AtomFmt::Bohr));
            }
        };
        if(type){
            // register/look-up type, then create bond
            this->bonds->bonds.push_back({at1, at2, getDistance(at1, at2, diff), diff,
                                            &*this->bonds->types.emplace(*type,
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
        bond.type = &*this->bonds->types.emplace(
                    type,
                    defaultColors[this->bonds->types.size()%5]
                    ).first;
    }

    // Modifier functions
    void modScale(AtomFmt tgt){
        if(tgt == this->at_fmt){ return; }
        this->evaluateCache();
        asFmt(tgt).evaluateCache();
        iterator source{this->atoms, this->pte, this->at_fmt, 0};
        iterator target{this->atoms, this->pte, tgt, 0};
        while(source.getIdx() != this->getNat()){
            target->coord = source->coord;
            ++target;
            ++source;
        }
        this->at_fmt = tgt;
    }
    void modShift(Vec shift, float fac=1.0f){
        shift *= fac;
        for(Atom& at:*this){
            at.coord += shift;
        }
    }
    void modRotate(float angle, Vec axis, Vec shift={0,0,0}){
        angle *= deg2rad;
        float c = std::cos(angle);
        float s = -std::sin(angle);
        float ic = 1.f-c;
        auto relative = this->at_fmt >= AtomFmt::Crystal;
        if(relative) axis = this->formatVec(axis, this->at_fmt, AtomFmt::Bohr);
        float len = Vec_length(axis);
        if(float_comp(len, 0.f)){
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
        if(!relative){
            for(Atom& at:*this){
                at.coord = (at.coord - shift) * rotMat + shift;
            }
        }else{
            const auto fwd = this->getFormatter(this->at_fmt, AtomFmt::Bohr);
            const auto bwd = this->getFormatter(AtomFmt::Bohr, this->at_fmt);
            shift = fwd(shift);
            for(Atom& at:*this){
                at.coord = bwd((fwd(at.coord) - shift) * rotMat + shift);
            }
        }
    }
    void modMirror(Vec ax1, Vec ax2, Vec shift={0,0,0}){
        auto relative = this->at_fmt >= AtomFmt::Crystal;
        if(relative){
            ax1 = this->formatVec(ax1, this->at_fmt, AtomFmt::Bohr);
            ax2 = this->formatVec(ax2, this->at_fmt, AtomFmt::Bohr);
            shift = this->formatVec(shift, this->at_fmt, AtomFmt::Bohr);
        }
        Vec normal = Vec_cross(ax1, ax2);
        normal /= Vec_length(normal);
        if(!relative){
            for(Atom& at:*this){
                at.coord -= 2*Vec_dot(at.coord-shift, normal)*normal;
            }
        }else{
            auto fwd = this->getFormatter(this->at_fmt, AtomFmt::Bohr);
            auto bwd = this->getFormatter(AtomFmt::Bohr, this->at_fmt);
            for(Atom& at:*this){
                at.coord -= bwd(2*Vec_dot(fwd(at.coord)-shift, normal)*normal);
            }
        }
    }

protected:
    StepMutable(std::shared_ptr<PeriodicTable> pte, AtomFmt fmt,
             std::shared_ptr<T> atoms, std::shared_ptr<BondList> bonds,
             std::shared_ptr<CellData> cell, std::shared_ptr<std::string> comment)
        : StepConst<T>{pte, fmt, atoms, bonds, cell, comment}
    {}
private:
    // Interface only:
    StepMutable(const StepMutable&) = default;
    StepMutable(StepMutable&&) = default;
    StepMutable& operator=(const StepMutable&) = default;
    StepMutable& operator=(StepMutable&&) = default;
};

}

#endif // STEPMUTABLE_H
