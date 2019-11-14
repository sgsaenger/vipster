#ifndef ORCA_PARAM_H
#define ORCA_PARAM_H

#include "../parameters.h"

namespace Vipster::IO{

struct OrcaParam final: BaseParam{
    using Header = std::vector<std::string>;
    OrcaParam(Header={});
    Header header;
    const struct Plugin* getFmt() const override;
    std::unique_ptr<BaseParam> copy() const override;
    void parseJson(const nlohmann::json&) override;
    nlohmann::json toJson() const override;
};

}

#endif // ORCA_PARAM_H
