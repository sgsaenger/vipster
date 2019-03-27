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
        return selection{this->pse, this->at_fmt, this, filter, this->cell, this->comment};
    }
    selection select(SelectionFilter filter)
    {
        return selection{this->pse, this->at_fmt, this, filter, this->cell, this->comment};
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
        auto tmp = StepMutable{this->pse, tgt,
                               this->atoms, this->bonds,
                               this->cell, this->comment};
        tmp.evaluateCache();
        return tmp;
    }

    // Atoms
    using StepConst<T>::operator[];
    Atom        operator[](size_t i) noexcept
    {
        return *iterator{this->atoms, this->at_fmt, i};
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
        return iterator{this->atoms, this->at_fmt, 0};
    }
    using StepConst<T>::end;
    iterator    end() noexcept
    {
        return iterator{this->atoms, this->at_fmt, this->getNat()};
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

    // Modifier functions
    void modScale(AtomFmt tgt){
        if(tgt == this->at_fmt){ return; }
        this->evaluateCache();
        asFmt(tgt).evaluateCache();
        iterator source{this->atoms, this->at_fmt, 0};
        iterator target{this->atoms, tgt, 0};
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
    StepMutable(std::shared_ptr<PseMap> pse, AtomFmt fmt,
             std::shared_ptr<T> atoms, std::shared_ptr<BondList> bonds,
             std::shared_ptr<CellData> cell, std::shared_ptr<std::string> comment)
        : StepConst<T>{pse, fmt, atoms, bonds, cell, comment}
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
