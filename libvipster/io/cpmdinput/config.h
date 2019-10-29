#ifndef CPI_CONF_H
#define CPI_CONF_H

#include "../configs.h"

namespace Vipster::IO{

struct CPConfig final: BaseConfig{
    enum class AtomFmt {Bohr, Angstrom, Crystal, Alat, Active};
    AtomFmt fmt;
    CPConfig(std::string="", AtomFmt=AtomFmt::Active);
    const struct Plugin* getFmt() const override;
    std::unique_ptr<BaseConfig> copy() const override;
    void parseJson(const nlohmann::json::iterator&) override;
    nlohmann::json toJson() const override;
};

void to_json(nlohmann::json& j, const CPConfig& p);
void from_json(const nlohmann::json& j, CPConfig& p);

}

#endif // CPI_CONF_H
