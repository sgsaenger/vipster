#ifndef CPI_PARAM_H
#define CPI_PARAM_H

#include "../parameters.h"

namespace Vipster::IO{

struct CPParam final: BaseParam{
    using Section = std::vector<std::string>;
    Section info;
    Section cpmd;
    Section system;
    Section pimd;
    Section path;
    Section ptddft;
    Section atoms;
    Section dft;
    Section prop;
    Section resp;
    Section linres;
    Section tddft;
    Section hardness;
    Section classic;
    Section exte;
    Section vdw;
    Section qmmm;
    std::string PPPrefix;
    std::string PPSuffix;
    std::string PPNonlocality;
    static const std::vector<std::pair<std::string, Section CPParam::*>> str2section;
    CPParam(Section={}, Section={}, Section={}, Section={},
            Section={}, Section={}, Section={}, Section={},
            Section={}, Section={}, Section={}, Section={},
            Section={}, Section={}, Section={}, Section={},
            Section={}, std::string="", std::string="", std::string="");
    const struct Plugin* getFmt() const override;
    std::unique_ptr<BaseParam> copy() const override;
    void parseJson(const nlohmann::json&) override;
    nlohmann::json toJson() const override;
};

}

#endif // CPI_PARAM_H
