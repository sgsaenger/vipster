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

inline std::ostream& operator<< (std::ostream& s, const KPoints &k)
{
    const KPoints::MPG &m = k.mpg;
    const KPoints::Discrete &d = k.discrete;
    switch (k.active) {
    case KPointFmt::Gamma:
        s << "Gamma-point only";
        break;
    case KPointFmt::MPG:
        s << "Monkhorst-Pack grid:\n"
          << " x: " << m.x
          << " y: " << m.y
          << " z: " << m.z
          << "\n sx: " << m.sx
          << " sy: " << m.sy
          << " sz: " << m.sz;
        break;
    case KPointFmt::Discrete:
        s << d.kpoints.size() << "discrete points:";
        if (d.properties & KPoints::band){
            s << "\n Band structure paths";
        }
        if (d.properties & KPoints::crystal){
            s << "\n Crystal coordinates";
        }
        for(const DiscreteKPoint &p:d.kpoints){
            s << "\n at: " << p.pos << ", Weight: " << p.weight;
        }
        break;
    }
    return s;
}
}

#endif // KPOINTS_H
