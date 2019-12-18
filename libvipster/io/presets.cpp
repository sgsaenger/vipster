#include "presets.h"
#include "../fileio.h"

using namespace Vipster;

IO::BasePreset::BasePreset(const struct Plugin* fmt,
                           CustomMap<std::string, PresetValue> &&values)
    : CustomMap{values}, fmt{fmt}
{}

const IO::Plugin* IO::BasePreset::getFmt() const
{
    return fmt;
}

void IO::to_json(nlohmann::json& j, const BasePreset& p)
{
    for(const auto &v: p){
        switch (v.second.index()) {
        case BasePreset::i_bool:
            j[v.first] = std::get<bool>(v.second);
            break;
        case BasePreset::i_uint:
            j[v.first] = std::get<uint>(v.second);
            break;
        default:
            break;
        }
    }
}

void IO::from_json(const nlohmann::json& j, BasePreset& p)
{
    for(auto &v: p){
        switch(v.second.index()) {
        case BasePreset::i_bool:
//            v.second = j.value(v.first, std::get<bool>(v.second));
        case BasePreset::i_uint:
//            v.second = j.value(v.first, std::get<uint>(v.second));
        default:
            break;
        }
    }
}
