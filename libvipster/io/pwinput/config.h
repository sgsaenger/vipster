#ifndef PWI_CONF_H
#define PWI_CONF_H

#include "../plugin.h"

namespace Vipster::IO{

struct PWConfig final: BaseConfig{
    enum class AtomFmt {Bohr, Angstrom, Crystal, Alat, Active};
    enum class CellFmt {Bohr, Angstrom, Active};
    AtomFmt atoms;
    CellFmt cell;
    PWConfig(std::string="", AtomFmt=AtomFmt::Active, CellFmt=CellFmt::Active);
    IOFmt getFmt() const override;
    std::unique_ptr<BaseConfig> copy() const override;
    void parseJson(const nlohmann::json::iterator&) override;
    nlohmann::json toJson() const override;
};

void to_json(nlohmann::json& j,const PWConfig& p);
void from_json(const nlohmann::json& j, PWConfig& p);

}

#endif // PWI_CONF_H
