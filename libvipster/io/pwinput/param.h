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
    std::string PPPrefix;
    std::string PPSuffix;
    static const std::map<std::string, Namelist PWParam::*> str2nl;
    PWParam(Namelist={}, Namelist={},
            Namelist={}, Namelist={}, Namelist={},
            std::string="", std::string="");
    IOFmt getFmt() const override;
    std::unique_ptr<BaseParam> copy() const override;
    void parseJson(const nlohmann::json&) override;
    nlohmann::json toJson() const override;
};

}

#endif // PWI_PARAM_H
