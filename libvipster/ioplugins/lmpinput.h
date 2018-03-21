#ifndef LMPINPUT_H
#define LMPINPUT_H

#include "../ioplugin.h"

namespace Vipster{
namespace IO{

extern const IOPlugin LmpInput;

enum class LmpAtomStyle{Angle, Atomic, Bond, Charge, Full, Molecular};

struct LmpConfig: BaseConfig{
    LmpAtomStyle style;
    bool angles, bonds, dihedrals, impropers;
};

}
}

#endif // LMPINPUT_H
