#ifndef CPI_CONF_H
#define CPI_CONF_H

#include "../plugin.h"

namespace Vipster::IO{

struct CPPreset final: BasePreset{
    enum class AtomFmt {Bohr, Angstrom, Crystal, Alat, Active};
    AtomFmt fmt;
    CPPreset(AtomFmt=AtomFmt::Active);
    const struct Plugin* getFmt() const override;
    std::unique_ptr<BasePreset> copy() const override;
    void parseJson(const nlohmann::json&) override;
    nlohmann::json toJson() const override;
};

}

#endif // CPI_CONF_H
