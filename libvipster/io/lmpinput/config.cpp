#include "config.h"

using namespace Vipster;

IO::LmpConfig::LmpConfig(std::string name, AtomStyle style,
                         bool bonds, bool angles, bool dihedrals, bool impropers)
    : BaseConfig{name},
      style{style},
      bonds{bonds}, angles{angles},
      dihedrals{dihedrals}, impropers{impropers}
{}

std::unique_ptr<BaseConfig> IO::LmpConfig::copy()
{
    return std::make_unique<IO::LmpConfig>(*this);
}
