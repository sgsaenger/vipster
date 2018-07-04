#include "param.h"

using namespace Vipster;

IO::PWParam::PWParam(std::string name, IO::PWParam::Namelist control, IO::PWParam::Namelist system,
                     IO::PWParam::Namelist electrons, IO::PWParam::Namelist ions,
                     IO::PWParam::Namelist cell)
    : BaseParam{name}, control{control}, system{system},
      electrons{electrons}, ions{ions}, cell{cell}
{}

std::unique_ptr<BaseParam> IO::PWParam::copy()
{
    return std::make_unique<IO::PWParam>(*this);
}

void Vipster::IO::to_json(nlohmann::json& j,const IO::PWParam& p)
{
    const std::map<std::string, IO::PWParam::Namelist IO::PWParam::*> nlmap = {
        {"&CONTROL", &IO::PWParam::control},
        {"&SYSTEM", &IO::PWParam::system},
        {"&ELECTRONS", &IO::PWParam::electrons},
        {"&IONS", &IO::PWParam::ions},
        {"&CELL", &IO::PWParam::cell},
    };
    j["name"] = p.name;
    for(auto& nl:nlmap){
        nlohmann::json& jnl = j.at(nl.first);
        for(auto kv:p.*nl.second){
            jnl[kv.first] = kv.second;
        }
    }
}

void Vipster::IO::from_json(const nlohmann::json& j, IO::PWParam& p)
{
    const std::map<std::string, IO::PWParam::Namelist IO::PWParam::*> nlmap = {
        {"&CONTROL", &IO::PWParam::control},
        {"&SYSTEM", &IO::PWParam::system},
        {"&ELECTRONS", &IO::PWParam::electrons},
        {"&IONS", &IO::PWParam::ions},
        {"&CELL", &IO::PWParam::cell},
    };
    p.name = j.at("name");
    for(auto& i:nlmap){
        IO::PWParam::Namelist& nl = p.*(i.second);
        for(auto it=j.at(i.first).begin(); it!=j.at(i.first).end(); ++it){
            nl[it.key()] = it.value();
        }
    }
}
