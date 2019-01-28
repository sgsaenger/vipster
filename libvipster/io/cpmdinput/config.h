#ifndef CPI_CONF_H
#define CPI_CONF_H

#include "../plugin.h"

namespace Vipster::IO{

struct CPConfig: BaseConfig{
    enum class AtomFmt {Bohr, Angstrom, Crystal, Alat, Current};
    AtomFmt fmt;
    CPConfig(std::string="", AtomFmt=AtomFmt::Current);
    std::unique_ptr<BaseConfig> copy() override;
    void parseJson(const nlohmann::json::iterator&) override;
    nlohmann::json toJson() override;
};

void to_json(nlohmann::json& j, const CPConfig& p);
void from_json(const nlohmann::json& j, CPConfig& p);

}

#endif // CPI_CONF_H
