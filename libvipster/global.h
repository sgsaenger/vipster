#ifndef GLOBAL_H
#define GLOBAL_H

namespace Vipster {
    constexpr float bohrrad = 0.52917721092;
    constexpr float invbohr = 1/bohrrad;

    enum class Fmt { Bohr = 1, Angstrom = 2, Crystal = 3, Alat = 4 };
}

#endif // GLOBAL_H
