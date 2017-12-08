#ifndef LMPINPUT_H
#define LMPINPUT_H

#include "ioplugin.h"

namespace Vipster{
namespace IO{
extern const IOPlugin LmpInput;
enum class LmpAtomStyle{Angle, Atomic, Bond, Charge, Full, Molecular};
struct LmpParam: BaseParam{
    LmpAtomStyle style;
    bool angles, bonds, dihedrals, impropers;
};
struct LmpData: Vipster::IO::BaseData{
    LmpParam data;
};
}
}

#endif // LMPINPUT_H
