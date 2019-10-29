#ifndef ORCA_PARAM_H
#define ORCA_PARAM_H

#include "../parameters.h"

namespace Vipster::IO{

struct OrcaParam final: BaseParam{
    using Header = std::vector<std::string>;
    OrcaParam(std::string="", Header={});
    Header header;
    const struct Plugin* getFmt() const override;
    std::unique_ptr<BaseParam> copy() const override;
    void parseJson(const nlohmann::json::iterator&) override;
    nlohmann::json toJson() const override;
};

void to_json(nlohmann::json& j, const OrcaParam& p);
void from_json(const nlohmann::json& j, OrcaParam& p);

}

#endif // ORCA_PARAM_H
