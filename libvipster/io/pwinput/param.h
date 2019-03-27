#ifndef PWI_PARAM_H
#define PWI_PARAM_H

#include "../plugin.h"

namespace Vipster::IO{

struct PWParam final: BaseParam{
    using Namelist = std::map<std::string, std::string>;
    Namelist control;
    Namelist system;
    Namelist electrons;
    Namelist ions;
    Namelist cell;
    static const std::map<std::string, Namelist PWParam::*> str2nl;
    PWParam(std::string="", Namelist={}, Namelist={},
            Namelist={}, Namelist={}, Namelist={});
    IOFmt getFmt() const override;
    std::unique_ptr<BaseParam> copy() const override;
    void parseJson(const nlohmann::json::iterator&) override;
    nlohmann::json toJson() const override;
};

void to_json(nlohmann::json& j,const PWParam& p);
void from_json(const nlohmann::json& j, PWParam& p);

}

#endif // PWI_PARAM_H
