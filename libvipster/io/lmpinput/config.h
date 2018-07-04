#ifndef LMPI_CONF_H
#define LMPI_CONF_H

#include "../plugin.h"

namespace Vipster{
namespace IO{

struct LmpConfig: BaseConfig{
    enum class AtomStyle{Angle, Atomic, Bond, Charge, Full, Molecular};
    AtomStyle style;
    bool bonds, angles, dihedrals, impropers;
    LmpConfig(std::string="", AtomStyle=AtomStyle::Atomic,
              bool=false, bool=false, bool=false, bool=false);
    std::unique_ptr<BaseConfig> copy() override;
};

const LmpConfig LmpConfigDefault{
    "default",
    LmpConfig::AtomStyle::Atomic,
    false, false, false, false
};

}
}

#endif // LMPI_CONF_H
