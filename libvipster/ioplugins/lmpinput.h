#ifndef LMPINPUT_H
#define LMPINPUT_H

#include "../ioplugin.h"

namespace Vipster{
namespace IO{

extern const IO::Plugin LmpInput;

enum class LmpAtomStyle{Angle, Atomic, Bond, Charge, Full, Molecular};

struct LmpConfig: BaseConfig{
    LmpAtomStyle style;
    bool angles, bonds, dihedrals, impropers;
    std::unique_ptr<BaseConfig> copy();
};

}
}

#endif // LMPINPUT_H
