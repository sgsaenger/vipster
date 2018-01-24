#ifndef ATOMREF_H
#define ATOMREF_H

#include "atom.h"

namespace Vipster {
    /*
     * Atom that references another Atom-like
     */
    class AtomRef: public Atom{
    public:
        AtomRef(const std::string *n, const Vec *co, const float *ch,
                const FixVec *f, const uint8_t *h, const bool *c_m, const bool *p_m);
        AtomRef(Atom& rhs);
        AtomRef& operator=(Atom& rhs);
    };
}

#endif // ATOMREF_H
