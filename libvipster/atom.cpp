#include <tuple>
#include <iostream>
#include "atom.h"

namespace Vipster{

Atom::Atom(const std::string *n, const Vec *co, const float *ch,
           const FixVec *f, const uint8_t *h, const bool *m)
    : name{n,m}, coord{co,m}, charge{ch,m}, fix{f,m}, hidden{h,m} {}

Atom& Atom::operator++()
{
    name.p_prop++;
    coord.p_prop++;
    charge.p_prop++;
    fix.p_prop++;
    hidden.p_prop++;
    return *this;
}

std::ostream& operator<<(std::ostream& s, const Atom::PropRef<std::string>& pr)
{
    s << static_cast<const std::string&>(pr);
    return s;
}
std::istream& operator>>(std::istream& s, Atom::PropRef<std::string>& pr)
{
    s >> static_cast<std::string&>(pr);
    return s;
}
bool operator==(const std::string& s, const Atom::PropRef<std::string>& pr)
{
    return s == static_cast<const std::string&>(pr);
}
bool operator==(const Atom::PropRef<std::string>& pr, const std::string& s)
{
    return s == static_cast<const std::string&>(pr);
}

//TODO: maybe ignore fix/hidden?
bool operator==(const Atom &a1, const Atom &a2) {
    return std::tie(a1.name, a1.coord, a1.charge, a1.fix, a1.hidden)
           ==
           std::tie(a2.name, a2.coord, a2.charge, a2.fix, a2.hidden);
}
bool operator!=(const Atom &a1, const Atom &a2) {
    return !(a1==a2);
}

}
