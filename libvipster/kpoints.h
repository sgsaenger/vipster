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
    KPointFmt active = KPointFmt::Gamma;
    struct MPG{
        int x,y,z;
        float sx,sy,sz;
    } mpg;
    struct Discrete{
        bool crystal,band;
        std::vector<DiscreteKPoint> kpoints;
    } discrete;
};

}

#endif // KPOINTS_H
