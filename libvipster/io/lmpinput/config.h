#ifndef LMPI_CONF_H
#define LMPI_CONF_H

#include "../plugin.h"

namespace Vipster::IO{

struct LmpConfig final: BaseConfig{
    enum class AtomStyle{Angle, Atomic, Bond, Charge, Full, Molecular};
    AtomStyle style;
    bool bonds, angles, dihedrals, impropers;
    const static std::map<AtomStyle, std::string> fmt2str;
    LmpConfig(std::string="", AtomStyle=AtomStyle::Atomic,
              bool=false, bool=false, bool=false, bool=false);
    IOFmt getFmt() const override;
    std::unique_ptr<BaseConfig> copy() const override;
    void parseJson(const nlohmann::json::iterator&) override;
    nlohmann::json toJson() const override;
};

void to_json(nlohmann::json& j,const LmpConfig& p);
void from_json(const nlohmann::json& j, LmpConfig& p);

}

#endif // LMPI_CONF_H
