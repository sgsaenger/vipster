#include "presets.h"
#include "fileio.h"
#include <nlohmann/json.hpp>

using namespace Vipster;

Preset::Preset(const struct Plugin* fmt, const BaseMap &values)
    : StaticMap{values}, fmt{fmt}
{}

const Plugin* Preset::getFmt() const
{
    return fmt;
}

namespace Vipster{
void to_json(nlohmann::json& j, const Preset& p)
{
    for(const auto &v: p){
        switch (v.second.first.index()) {
        case Preset::i_bool:
            j[v.first] = std::get<bool>(v.second.first);
            break;
        case Preset::i_enum:
            j[v.first] = static_cast<const std::string&>(std::get<NamedEnum>(v.second.first));
            break;
        default:
            break;
        }
    }
}

void from_json(const nlohmann::json& j, Preset& p)
{
    for(auto &v: p){
        switch(v.second.first.index()) {
        case Preset::i_bool:
            v.second.first = j.value(v.first, std::get<bool>(v.second.first));
            break;
        case Preset::i_enum:
        {
            auto &val = std::get<NamedEnum>(v.second.first);
            if(auto pos = j.find(v.first); pos != j.end()){
                try{
                    val = pos->get<std::string>();
                }catch(IOError&){
                    // ignore
                }
            }
        }
            break;
        default:
            break;
        }
    }
}
}
