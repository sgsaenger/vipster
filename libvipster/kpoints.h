#ifndef KPOINTS_H
#define KPOINTS_H

#include "vec.h"
#include <vector>

namespace Vipster {

struct KPoints{
    enum class Fmt{Gamma, MPG, Discrete};
    Fmt active = Fmt::Gamma;
    struct MPG{
        int x{1},y{1},z{1};
        float sx{},sy{},sz{};
    } mpg;
    struct Discrete{
        enum Properties{none=0x0,crystal=0x1,band=0x2,contour=0x4};
        struct Point{
            Vec pos;
            float weight;
        };
        uint8_t properties{};
        std::vector<Point> kpoints{};
    } discrete;
};
}

#endif // KPOINTS_H
