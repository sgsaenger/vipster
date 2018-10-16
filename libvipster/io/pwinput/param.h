#ifndef PWI_PARAM_H
#define PWI_PARAM_H

#include "../plugin.h"

namespace Vipster{
namespace IO{

struct PWParam: BaseParam{
    using Namelist = std::map<std::string, std::string>;
    Namelist control;
    Namelist system;
    Namelist electrons;
    Namelist ions;
    Namelist cell;
    static const std::map<std::string, Namelist PWParam::*> str2nl;
    PWParam(std::string="", Namelist={}, Namelist={},
            Namelist={}, Namelist={}, Namelist={});
    std::unique_ptr<BaseParam> copy() override;
    void parseJson(const nlohmann::json::iterator&) override;
    nlohmann::json toJson() override;
};

void to_json(nlohmann::json& j,const PWParam& p);
void from_json(const nlohmann::json& j, PWParam& p);

}
}

#endif // PWI_PARAM_H
