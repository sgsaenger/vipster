#include "step.h"
#include <algorithm>

using namespace Vipster;

Step::Step(std::shared_ptr<PseMap> pse, AtomFmt fmt)
    : pse{pse}, fmt{fmt} {}

AtomFmt Step::getFmt() const noexcept
{
    return fmt;
}

size_t Step::getNat() const noexcept
{
    return at_coord->size();
}

Step::iterator Step::begin()
{
    return Step::iterator(this, 0);
}

const Step::iterator Step::begin() const
{
    return Step::iterator(this, 0);
}

const Step::iterator Step::cbegin() const
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

const Step::iterator Step::cend() const
{
    return Step::iterator(this, getNat());
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

std::function<Vec(Vec)> Step::getFormatter(AtomFmt source, AtomFmt target) const noexcept
{
    //TODO: getInvDim?
    switch(source) {
    case AtomFmt::Bohr:
        switch(target){
        case AtomFmt::Angstrom:
            return [](Vec v){v*Vipster::bohrrad;return v;};
        case AtomFmt::Alat:
            return [this](Vec v){v/getCellDim(CdmFmt::Bohr);return v;};
        case AtomFmt::Crystal:
            return [this](Vec v){v*getInvVec()/getCellDim(CdmFmt::Bohr);return v;};
        default:
            break;
        }
        break;
    case AtomFmt::Angstrom:
        switch(target){
        case AtomFmt::Bohr:
            return [](Vec v){v*Vipster::invbohr;return v;};
        case AtomFmt::Alat:
            return [this](Vec v){v/getCellDim(CdmFmt::Angstrom);return v;};
        case AtomFmt::Crystal:
            return [this](Vec v){v*getInvVec()/getCellDim(CdmFmt::Angstrom);return v;};
        default:
            break;
        }
        break;
    case AtomFmt::Alat:
        switch(target){
        case AtomFmt::Angstrom:
            return [this](Vec v){v*getCellDim(CdmFmt::Angstrom);return v;};
        case AtomFmt::Bohr:
            return [this](Vec v){v*getCellDim(CdmFmt::Bohr);return v;};
        case AtomFmt::Crystal:
            return [this](Vec v){v*getInvVec();return v;};
        default:
            break;
        }
        break;
    case AtomFmt::Crystal:
        switch(target){
        case AtomFmt::Angstrom:
            return [this](Vec v){v*getCellVec()*getCellDim(CdmFmt::Angstrom); return v;};
        case AtomFmt::Alat:
            return [this](Vec v){v*getCellVec(); return v;};
        case AtomFmt::Bohr:
            return [this](Vec v){v*getCellVec()*getCellDim(CdmFmt::Bohr); return v;};
        default:
            break;
        }
    }
    return [](Vec v){return v;};
}

Vec Step::formatVec(Vec in, AtomFmt source, AtomFmt target) const
{
    return getFormatter(source, target)(in);
}

std::vector<Vec> Step::formatAll(std::vector<Vec> in, AtomFmt source, AtomFmt target) const
{
    if ((source == target) || (in.size() == 0)) return in;
    auto op = getFormatter(source, target);
    std::transform(in.begin(), in.end(), in.begin(), op);
    return in;
}

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
