#include "atomref.h"

using namespace Vipster;

AtomRef::AtomRef(const std::string *n, const Vec *co, const float *ch,
                const FixVec *f, const uint8_t *h, const bool *c_m, const bool *p_m)
    : Atom{n, co, ch, f, h, c_m, p_m} {}

// like a real reference, constructing makes it point to the origin of rhs
AtomRef::AtomRef(Atom& rhs)
    : Atom{rhs} {
}

// also like a real reference, assigning changes the origin
AtomRef& AtomRef::operator=(Atom& rhs) {
    name = rhs.name;
    coord = rhs.coord;
    charge = rhs.charge;
    fix = rhs.fix;
    hidden = rhs.hidden;
    return *this;
}
