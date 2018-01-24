#include "step.h"
#include <algorithm>

using namespace Vipster;

Step::Step(std::shared_ptr<PseMap> pse, AtomFmt fmt)
    : pse{pse}, at_fmt{fmt} {}

AtomFmt Step::getFmt() const noexcept
{
    return at_fmt;
}

size_t Step::getNat() const noexcept
{
    evaluateCache();
    return at_coord->size();
}

Step::iterator Step::begin()
{
    evaluateCache();
    return Step::iterator(this, 0);
}

const Step::iterator Step::begin() const
{
    evaluateCache();
    return Step::iterator(this, 0);
}

const Step::iterator Step::cbegin() const
{
    evaluateCache();
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
    evaluateCache();
    std::set<std::string> types;
    for(const auto& at:*this) {
        types.insert(at.name);
    }
    return types;
}

size_t Step::getNtyp() const noexcept
{
    return getTypes().size();
}

Step::iterator::iterator(const Step *step, size_t idx)
    : AtomRef{(*step)[idx]}, step{const_cast<Step*>(step)}, idx{idx} {}

Step::iterator& Step::iterator::operator++()
{
    ++idx;
    Atom::operator++();
    return *this;
}

Step::iterator Step::iterator::operator++(int)
{
    Step::iterator tmp{*this};
    ++idx;
    Atom::operator++();
    return tmp;
}

AtomRef& Step::iterator::operator*()
{
    return static_cast<AtomRef&>(*this);
}

const AtomRef& Step::iterator::operator*() const
{
    return static_cast<const AtomRef&>(*this);
}

bool Step::iterator::operator==(const Step::iterator& rhs) const
{
    return ((step == rhs.step) || (idx == rhs.idx));
}

bool Step::iterator::operator!=(const Step::iterator& rhs) const
{
    return ((step != rhs.step) || (idx != rhs.idx));
}

std::function<Vec (const Vec &)> Step::getFormatter(AtomFmt source, AtomFmt target) const noexcept
{
    //TODO: getInvDim?
    //TODO: save/catch variables?
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
