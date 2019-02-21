#include "param.h"

using namespace Vipster;

IO::PWParam::PWParam(IO::PWParam::Namelist control,
                     IO::PWParam::Namelist system, IO::PWParam::Namelist electrons,
                     IO::PWParam::Namelist ions, IO::PWParam::Namelist cell,
                     std::string PPPrefix, std::string PPSuffix)
    : control{control}, system{system},
      electrons{electrons}, ions{ions}, cell{cell},
      PPPrefix{PPPrefix}, PPSuffix{PPSuffix}
{}

IOFmt IO::PWParam::getFmt() const
{
    return IOFmt::PWI;
}

std::unique_ptr<IO::BaseParam> IO::PWParam::copy() const
{
    return std::make_unique<IO::PWParam>(*this);
}

const std::map<std::string, IO::PWParam::Namelist IO::PWParam::*> IO::PWParam::str2nl = {
    {"&CONTROL", &IO::PWParam::control},
    {"&SYSTEM", &IO::PWParam::system},
    {"&ELECTRONS", &IO::PWParam::electrons},
    {"&IONS", &IO::PWParam::ions},
    {"&CELL", &IO::PWParam::cell},
};

void IO::PWParam::parseJson(const nlohmann::json &j)
{
    for(auto& nl: PWParam::str2nl){
        this->*nl.second = j.value(nl.first, PWParam::Namelist{});
    }
    PPPrefix = j.value("PPPrefix", "");
    PPSuffix = j.value("PPSuffix", "");
}

nlohmann::json IO::PWParam::toJson() const
{
    nlohmann::json j;
    for(auto& nl: PWParam::str2nl){
        j[nl.first] = this->*nl.second;
    }
    j["PPPrefix"] = PPPrefix;
    j["PPSuffix"] = PPSuffix;
    return j;
}
