#include "parameters.h"
#include "fileio.h"

#include <nlohmann/json.hpp>

using namespace Vipster;

Parameter::Parameter(const struct Plugin* fmt, const BaseMap &values)
    : StaticMap{values}, fmt{fmt}
{}

const Plugin* Parameter::getFmt() const
{
    return fmt;
}

namespace Vipster{
void to_json(nlohmann::json& j, const Parameter& p)
{
    for(const auto &v: p){
        switch (v.second.first.index()) {
        case Parameter::i_str:
            j[v.first] = std::get<Parameter::i_str>(v.second.first);
            break;
        case Parameter::i_strvec:
            j[v.first] = std::get<Parameter::i_strvec>(v.second.first);
            break;
        case Parameter::i_strmap:
            j[v.first] = std::get<Parameter::i_strmap>(v.second.first);
            break;
        default:
            break;
        }
    }
}

void from_json(const nlohmann::json& j, Parameter& p)
{
    for(auto &v: p){
        switch (v.second.first.index()) {
        case Parameter::i_str:
            v.second.first = j.value(v.first,
                 std::get<Parameter::i_str>(v.second.first));
            break;
        case Parameter::i_strvec:
            v.second.first = j.value(v.first,
                 std::get<Parameter::i_strvec>(v.second.first));
            break;
        case Parameter::i_strmap:
            v.second.first = j.value(v.first,
                 std::get<Parameter::i_strmap>(v.second.first));
            break;
        default:
            break;
        }
    }
}
}
