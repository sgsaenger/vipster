#ifndef PWINPUT_H
#define PWINPUT_H

#include "../ioplugin.h"

namespace Vipster{
namespace IO{

extern const IO::Plugin PWInput;

using PWNamelist = std::map<std::string, std::string>;

struct PWParam: BaseParam{
    PWNamelist control;
    PWNamelist system;
    PWNamelist electrons;
    PWNamelist ions;
    PWNamelist cell;
    PWParam(std::string="", PWNamelist={}, PWNamelist={},
            PWNamelist={}, PWNamelist={}, PWNamelist={});
    std::unique_ptr<BaseParam> copy() override;
};

const PWParam PWParamDefault{
    "default",
    {}, {}, {}, {}, {}
};

void to_json(nlohmann::json& j,const PWParam& p);
void from_json(const nlohmann::json& j, PWParam& p);
//TODO: json methods for config

struct PWConfig: BaseConfig{
    enum class AtomFmt {Bohr, Angstrom, Crystal, Alat, Current};
    enum class CellFmt {Bohr, Angstrom, Current};
    AtomFmt atoms;
    CellFmt cell;
    PWConfig(std::string="", AtomFmt=AtomFmt::Current, CellFmt=CellFmt::Current);
    std::unique_ptr<BaseConfig> copy() override;
};

const PWConfig PWConfigDefault{
    "default",
    PWConfig::AtomFmt::Bohr,
    PWConfig::CellFmt::Bohr
};

}
}


#endif // PWINPUT_H
