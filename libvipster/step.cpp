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
    idx++;
    Atom::operator ++();
    return *this;
}

Step::iterator Step::iterator::operator++(int)
{
    Step::iterator tmp{*this};
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

std::function<Vec(Vec)> Step::getFormatter(AtomFmt source, AtomFmt target) const noexcept
{
    //TODO: getInvDim?
    //TODO: save/catch variables?
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
