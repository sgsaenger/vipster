#ifndef POSCAR_CONFIG_H
#define POSCAR_CONFIG_H

#include "../plugin.h"

namespace Vipster::IO{

struct PoscarConfig final: BaseConfig{
    bool selective;
    bool cartesian;
    PoscarConfig(std::string="", bool selective=true, bool cartesian=false);
    IOFmt getFmt() const override;
    std::unique_ptr<BaseConfig> copy() const override;
    void parseJson(const nlohmann::json::iterator &) override;
    nlohmann::json toJson() const override;
};

void to_json(nlohmann::json& j, const PoscarConfig& c);
void from_json(const nlohmann::json& j, PoscarConfig& c);

}

#endif // POSCAR_CONFIG_H
