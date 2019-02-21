#ifndef LMPI_CONF_H
#define LMPI_CONF_H

#include "../plugin.h"

namespace Vipster::IO{

struct LmpPreset final: BasePreset{
    enum class AtomStyle{Angle, Atomic, Bond, Charge, Full, Molecular};
    AtomStyle style;
    bool bonds, angles, dihedrals, impropers;
    const static std::map<AtomStyle, std::string> fmt2str;
    LmpPreset(AtomStyle=AtomStyle::Atomic,
              bool=false, bool=false, bool=false, bool=false);
    IOFmt getFmt() const override;
    std::unique_ptr<BasePreset> copy() const override;
    void parseJson(const nlohmann::json&) override;
    nlohmann::json toJson() const override;
};

}

#endif // LMPI_CONF_H
