#ifndef PWI_CONF_H
#define PWI_CONF_H

#include "../plugin.h"

namespace Vipster{
namespace IO{

//TODO: json methods for config

struct PWConfig: BaseConfig{
    enum class AtomFmt {Bohr, Angstrom, Crystal, Alat, Current};
    enum class CellFmt {Bohr, Angstrom, Current};
    AtomFmt atoms;
    CellFmt cell;
    PWConfig(std::string="", AtomFmt=AtomFmt::Current, CellFmt=CellFmt::Current);
    std::unique_ptr<BaseConfig> copy() override;
};

const PWConfig PWConfigDefault{
    "default",
    PWConfig::AtomFmt::Current,
    PWConfig::CellFmt::Current
};

}
}

#endif // PWI_CONF_H
