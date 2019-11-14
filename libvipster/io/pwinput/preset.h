#ifndef PWI_CONF_H
#define PWI_CONF_H

#include "../presets.h"

namespace Vipster::IO{

struct PWPreset final: BasePreset{
    enum class AtomFmt {Bohr, Angstrom, Crystal, Alat, Active};
    enum class CellFmt {Bohr, Angstrom, Active};
    AtomFmt atoms;
    CellFmt cell;
    PWPreset(AtomFmt=AtomFmt::Active, CellFmt=CellFmt::Active);
    const struct Plugin* getFmt() const override;
    std::unique_ptr<BasePreset> copy() const override;
    void parseJson(const nlohmann::json&) override;
    nlohmann::json toJson() const override;
};

}

#endif // PWI_CONF_H
