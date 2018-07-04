#ifndef PWINPUT_H
#define PWINPUT_H

#include "../ioplugin.h"

namespace Vipster{
namespace IO{

extern const IO::Plugin PWInput;

struct PWParam: BaseParam{
    using Namelist = std::map<std::string, std::string>;
    Namelist control;
    Namelist system;
    Namelist electrons;
    Namelist ions;
    Namelist cell;
    PWParam(std::string="", Namelist={}, Namelist={},
            Namelist={}, Namelist={}, Namelist={});
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
    PWConfig::AtomFmt::Current,
    PWConfig::CellFmt::Current
};

}
}


#endif // PWINPUT_H
