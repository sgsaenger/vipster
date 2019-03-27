#include "param.h"

using namespace Vipster;

IO::PWParam::PWParam(std::string name, IO::PWParam::Namelist control,
                     IO::PWParam::Namelist system, IO::PWParam::Namelist electrons,
                     IO::PWParam::Namelist ions, IO::PWParam::Namelist cell)
    : BaseParam{name}, control{control}, system{system},
      electrons{electrons}, ions{ions}, cell{cell}
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

void IO::to_json(nlohmann::json& j,const PWParam& p)
{
    for(auto& nl: PWParam::str2nl){
        j[nl.first] = p.*nl.second;
    }
}

void IO::from_json(const nlohmann::json& j, PWParam& p)
{
    for(auto& i: PWParam::str2nl){
        p.*i.second = j.value(i.first, PWParam::Namelist{});
    }
}

void IO::PWParam::parseJson(const nlohmann::json::iterator& it)
{
    name = it.key();
    from_json(it.value(), *this);
}

nlohmann::json IO::PWParam::toJson() const
{
    nlohmann::json j;
    to_json(j, *this);
    return j;
}
