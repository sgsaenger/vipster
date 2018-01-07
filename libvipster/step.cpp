#include "step.h"

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
