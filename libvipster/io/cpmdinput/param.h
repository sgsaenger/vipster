#ifndef CPI_PARAM_H
#define CPI_PARAM_H

#include "../plugin.h"

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
    static const std::vector<std::pair<std::string, Section CPParam::*>> str2section;
    CPParam(std::string="", Section={}, Section={}, Section={}, Section={},
            Section={}, Section={}, Section={}, Section={}, Section={}, Section={},
            Section={}, Section={}, Section={}, Section={}, Section={}, Section={});
    IOFmt getFmt() const override;
    std::unique_ptr<BaseParam> copy() const override;
    void parseJson(const nlohmann::json::iterator&) override;
    nlohmann::json toJson() const override;
};

void to_json(nlohmann::json& j,const CPParam& p);
void from_json(const nlohmann::json& j, CPParam& p);

}

#endif // CPI_PARAM_H
