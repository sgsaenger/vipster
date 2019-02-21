#include "presets.h"
#include "../io.h"

using namespace Vipster;

IO::Presets Vipster::presets = [](){
    IO::Presets tmp;
    for(const auto& p: IOPlugins){
        if(p.second->arguments & IO::Plugin::Args::Preset){
            tmp[p.first]["default"] = p.second->makePreset();
        }
    }
    return tmp;
}();

void IO::to_json(nlohmann::json& j, const BasePreset& p)
{
    j = p.toJson();
}

void IO::from_json(const nlohmann::json& j, BasePreset& p)
{
    p.parseJson(j);
}
