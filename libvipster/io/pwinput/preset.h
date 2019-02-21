#ifndef PWI_CONF_H
#define PWI_CONF_H

#include "../plugin.h"

namespace Vipster::IO{

struct PWPreset final: BasePreset{
    enum class AtomFmt {Bohr, Angstrom, Crystal, Alat, Current};
    enum class CellFmt {Bohr, Angstrom, Current};
    AtomFmt atoms;
    CellFmt cell;
    PWPreset(AtomFmt=AtomFmt::Current, CellFmt=CellFmt::Current);
    IOFmt getFmt() const override;
    std::unique_ptr<BasePreset> copy() const override;
    void parseJson(const nlohmann::json&) override;
    nlohmann::json toJson() const override;
};

}

#endif // PWI_CONF_H
