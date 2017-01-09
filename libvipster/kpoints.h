#ifndef KPOINTS_H
#define KPOINTS_H

#include <vec.h>
#include <vector>

namespace Vipster {
enum class KPointFmt{Gamma, MPG, Discrete};
struct DiscreteKPoint{
    Vec pos;
    float weight;
};

struct KPoints{
    enum DiscreteProperties{none=0x0,crystal=0x1,band=0x2};
    KPointFmt active = KPointFmt::Gamma;
    struct MPG{
        int x,y,z;
        float sx,sy,sz;
    } mpg;
    struct Discrete{
        DiscreteProperties properties;
        std::vector<DiscreteKPoint> kpoints;
    } discrete;
};

}

#endif // KPOINTS_H
