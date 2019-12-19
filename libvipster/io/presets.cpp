#include "presets.h"
#include "../fileio.h"

using namespace Vipster;

IO::Preset::Preset(const struct Plugin* fmt,
                           CustomMap<std::string, PresetValue> &&values)
    : CustomMap{values}, fmt{fmt}
{}

const IO::Plugin* IO::Preset::getFmt() const
{
    return fmt;
}

void IO::to_json(nlohmann::json& j, const Preset& p)
{
    for(const auto &v: p){
        switch (v.second.index()) {
        case Preset::i_bool:
            j[v.first] = std::get<Preset::i_bool>(v.second);
            break;
        case Preset::i_uint:
            j[v.first] = std::get<Preset::i_uint>(v.second);
            break;
        default:
            break;
        }
    }
}

void IO::from_json(const nlohmann::json& j, Preset& p)
{
    for(auto &v: p){
        switch(v.second.index()) {
        case Preset::i_bool:
            v.second = j.value(v.first, std::get<Preset::i_bool>(v.second));
            break;
        case Preset::i_uint:
            v.second = j.value(v.first, std::get<Preset::i_uint>(v.second));
            break;
        default:
            break;
        }
    }
}
