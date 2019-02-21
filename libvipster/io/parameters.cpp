#include "parameters.h"
#include "../io.h"

using namespace Vipster;

IO::Parameters Vipster::params = [](){
    IO::Parameters tmp;
    for(const auto& p: IOPlugins){
        if(p.second->arguments & IO::Plugin::Args::Param){
            tmp[p.first]["default"] = p.second->makeParam();
        }
    }
    return tmp;
}();

void IO::to_json(nlohmann::json& j, const BaseParam& p)
{
    j = p.toJson();
}

void IO::from_json(const nlohmann::json& j, BaseParam& p)
{
    p.parseJson(j);
}
