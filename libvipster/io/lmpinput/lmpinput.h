#ifndef LMPINPUT_H
#define LMPINPUT_H

#include "../plugin.h"

namespace Vipster{
namespace IO{

extern const IO::Plugin LmpInput;

enum class LmpAtomStyle{Angle, Atomic, Bond, Charge, Full, Molecular};

struct LmpConfig: BaseConfig{
    LmpAtomStyle style;
    bool bonds, angles, dihedrals, impropers;
    LmpConfig(std::string="", LmpAtomStyle=LmpAtomStyle::Atomic,
              bool=false, bool=false, bool=false, bool=false);
    std::unique_ptr<BaseConfig> copy() override;
};

const LmpConfig LmpConfigDefault{
    "default",
    LmpAtomStyle::Atomic,
    false, false, false, false
};

}
}

#endif // LMPINPUT_H
