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
    virtual ~StepMutable()=default;

    using iterator = typename T::iterator;
    using constIterator = typename T::constIterator;

    using StepConst<T>::asFmt;
    StepMutable&    asFmt(AtomFmt)
    {
        //TODO: BROKEN, IMPLEMENT!
        return *this;
    }

    // Comment
    void                setComment(const std::string& s)
    {
        *this->comment = s;
    }

    // Atoms
    using StepConst<T>::operator[];
    Atom        operator[](size_t i) noexcept
    {
        return *iterator{this->atoms, this->at_fmt, i};
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

    // Modifier functions
    void modScale(AtomFmt tgt){
        if(tgt == this->at_fmt){ return; }
        this->evaluateCache();
        asFmt(tgt).evaluateCache();
        typename T::iterator source{this->atoms, this->at_fmt, 0};
        typename T::iterator target{this->atoms, tgt, 0};
        while(source.getIdx() != this->getNat()){
            target->coord = source->coord;
            ++target;
            ++source;
        }
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
        for(Atom& at:*this){
            at.coord = (at.coord - shift) * rotMat + shift;
        }
    }
    void modMirror(Vec ax1, Vec ax2, Vec shift={0,0,0}){
        Vec normal = Vec_cross(ax1, ax2);
        normal /= Vec_length(normal);
        for(Atom& at:*this){
            at.coord -= 2*Vec_dot(at.coord-shift, normal)*normal;
        }
    }
    void modWrap(){
        for(Atom& at:asFmt(AtomFmt::Crystal)){
            at.coord[0] -= std::floor(at.coord[0]);
            at.coord[1] -= std::floor(at.coord[1]);
            at.coord[2] -= std::floor(at.coord[2]);
        }
    }
    void modCrop(){
        std::vector<size_t> toRemove;
        toRemove.reserve(this->getNat());
        for(auto it=asFmt(AtomFmt::Crystal).begin(); it!=asFmt(AtomFmt::Crystal).end(); ++it){
            if((it->coord[0]>=1) | (it->coord[0]<0) |
               (it->coord[1]>=1) | (it->coord[1]<0) |
               (it->coord[2]>=1) | (it->coord[2]<0)){
                toRemove.push_back(it.getIdx());
            }
        }
        for(auto it=toRemove.rbegin(); it!=toRemove.rend(); ++it){
            static_cast<T*>(this)->delAtom(*it);
        }
    }
    void modMultiply(size_t x, size_t y, size_t z){
        auto fac = x*y*z;
        if(fac == 0){
            throw Error("Cannot eradicate atoms via modMultiply");
        }else if(fac == 1){
            return;
        }
        auto& handle = asFmt(AtomFmt::Crystal);
        auto cell = this->getCellVec();
        auto multiply = [&](uint8_t dir, uint8_t mult){
            auto atoms = handle.getAtoms();
            auto oldNat = handle.getNat();
            cell[dir] *= mult;
            for(uint8_t i=1; i<mult; ++i){
                handle.newAtoms(atoms);
                auto refIt = handle.begin();
                for(auto it=refIt+i*oldNat; it!=refIt+(i+1)*oldNat; ++it){
                    it->coord[dir] += i;
                }
            }
        };
        if(x>1){
            multiply(0, x);
        }
        if(y>1){
            multiply(1, y);
        }
        if(z>1){
            multiply(2, z);
        }
        static_cast<T*>(this)->setCellVec(cell);
    }
    void modAlign(uint8_t step_dir, uint8_t target_dir){
        auto target = Vec{};
        target.at(target_dir) = 1;
        auto source = this->getCellVec().at(step_dir);
        source /= Vec_length(source);
        if(target == source){
            return;
        }
        auto axis = Vec_cross(source, target);
        axis /= Vec_length(axis);
        auto cos = Vec_dot(source, target);
        auto icos = 1-cos;
        auto sin = -std::sqrt(1-cos*cos);
        Mat rotMat = {Vec{icos * axis[0] * axis[0] + cos,
                          icos * axis[0] * axis[1] - sin * axis[2],
                          icos * axis[0] * axis[2] + sin * axis[1]},
                      Vec{icos * axis[1] * axis[0] + sin * axis[2],
                          icos * axis[1] * axis[1] + cos,
                          icos * axis[1] * axis[2] - sin * axis[0]},
                      Vec{icos * axis[2] * axis[0] - sin * axis[1],
                          icos * axis[2] * axis[1] + sin * axis[0],
                          icos * axis[2] * axis[2] + cos}};
        Mat oldCell = this->getCellVec();
        Mat newCell = oldCell*rotMat;
        static_cast<T*>(this)->setCellVec(newCell, true);
    }
    void modReshape(Mat newMat, float newCdm, CdmFmt cdmFmt){
        auto oldCdm = this->getCellDim(cdmFmt);
        auto oldMat = this->getCellVec();
        if((newMat == oldMat) && (float_comp(newCdm, oldCdm))){
            return;
        }
        modWrap();
        auto& handle = *static_cast<T*>(this);
        size_t fac;
        if(newMat == oldMat){
            // only changing cdm
            fac = std::ceil(newCdm/oldCdm);
        }else{
            // change vectors or both
            auto getExtent = [](const Mat& m){
                return Vec{m[0][0] + m[1][0] + m[2][0],
                           m[0][1] + m[1][1] + m[2][1],
                           m[0][2] + m[1][2] + m[2][2]
                           };
            };
            auto compExtLt = [](const Vec& v1, const Vec& v2){
                return (v1[0] < v2[0]) || (v1[1] < v2[1]) || (v1[2] < v2[2]);
            };
            Vec newExtent = getExtent(newMat*newCdm);
            Vec oldExtent = getExtent(oldMat*oldCdm);
            fac = 1;
            while(compExtLt(oldExtent*fac, newExtent)){
                fac += 1;
            }
        }
        modMultiply(fac, fac, fac);
        handle.setCellVec(newMat);
        handle.setCellDim(newCdm, cdmFmt);
        modCrop();
    }

protected:
    StepMutable(std::shared_ptr<PseMap> pse, AtomFmt fmt,
             std::shared_ptr<T> atoms, std::shared_ptr<BondList> bonds,
             std::shared_ptr<CellData> cell, std::shared_ptr<std::string> comment)
        : StepConst<T>{pse, fmt, atoms, bonds, cell, comment}
    {}
};

}

#endif // STEPMUTABLE_H
