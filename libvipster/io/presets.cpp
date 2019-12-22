#include "presets.h"
#include "../fileio.h"

using namespace Vipster;

IO::Preset::Preset(const struct Plugin* fmt, BaseMap &&values)
    : StaticMap{values}, fmt{fmt}
{}

const IO::Plugin* IO::Preset::getFmt() const
{
    return fmt;
}

void IO::to_json(nlohmann::json& j, const Preset& p)
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

void IO::from_json(const nlohmann::json& j, Preset& p)
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
                }catch(IO::Error){
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
