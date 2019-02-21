#ifndef CPI_CONF_H
#define CPI_CONF_H

#include "../plugin.h"

namespace Vipster::IO{

struct CPPreset final: BasePreset{
    enum class AtomFmt {Bohr, Angstrom, Crystal, Alat, Current};
    AtomFmt fmt;
    CPPreset(AtomFmt=AtomFmt::Current);
    IOFmt getFmt() const override;
    std::unique_ptr<BasePreset> copy() const override;
    void parseJson(const nlohmann::json&) override;
    nlohmann::json toJson() const override;
};

}

#endif // CPI_CONF_H
