#ifndef POSCAR_CONFIG_H
#define POSCAR_CONFIG_H

#include "../presets.h"

namespace Vipster::IO{

struct PoscarPreset final: BasePreset{
    bool selective;
    bool cartesian;
    PoscarPreset(bool selective=true, bool cartesian=false);
    const struct Plugin* getFmt() const override;
    std::unique_ptr<BasePreset> copy() const override;
    void parseJson(const nlohmann::json &) override;
    nlohmann::json toJson() const override;
};

}

#endif // POSCAR_CONFIG_H
