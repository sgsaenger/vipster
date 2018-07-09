#ifndef CPI_PARAM_H
#define CPI_PARAM_H

#include "../plugin.h"

namespace Vipster {
namespace IO {

struct CPParam: BaseParam{
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
    static const std::map<std::string, Section CPParam::*> str2section;
    CPParam(std::string="", Section={}, Section={}, Section={}, Section={},
            Section={}, Section={}, Section={}, Section={}, Section={}, Section={},
            Section={}, Section={}, Section={}, Section={}, Section={}, Section={});
    std::unique_ptr<BaseParam> copy() override;
    void parseJson(const nlohmann::json&) override;
    nlohmann::json toJson() override;
};

void to_json(nlohmann::json& j,const CPParam& p);
void from_json(const nlohmann::json& j, CPParam& p);

}
}

#endif // CPI_PARAM_H
