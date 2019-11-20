#ifndef CELL_H
#define CELL_H

#include "vec.h"
#include "global.h"
#include <vector>

namespace Vipster {
    enum class CdmFmt { Bohr, Angstrom };

    struct CellData {
        bool    enabled{false};
        double  dimBohr{1};
        double  dimAngstrom{bohrrad};
        Mat     cellvec{{{{1,0,0}},{{0,1,0}},{{0,0,1}}}};
        Mat     invvec{{{{1,0,0}},{{0,1,0}},{{0,0,1}}}};
    };
}

#endif // CELL_H
