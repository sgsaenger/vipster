#ifndef XYZ_CONF_H
#define XYZ_CONF_H

#include "../configs.h"

namespace Vipster::IO{

struct XYZConfig final: BaseConfig{
    enum class Data{None, Charge, Forces};
    enum class Mode{Step, Trajec, Cell};
    Mode filemode;
    Data atomdata;
    XYZConfig(std::string="", Mode=Mode::Step, Data=Data::None);
    const struct Plugin* getFmt() const override;
    std::unique_ptr<BaseConfig> copy() const override;
    void parseJson(const nlohmann::json::iterator&) override;
    nlohmann::json toJson() const override;
};

void to_json(nlohmann::json& j, const XYZConfig& p);
void from_json(const nlohmann::json& j, XYZConfig& p);

}

#endif // XYZ_CONF_H
