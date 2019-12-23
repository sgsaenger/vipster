#include "parameters.h"
#include "../fileio.h"

using namespace Vipster;

IO::Parameter::Parameter(const struct Plugin* fmt, const BaseMap &values)
    : StaticMap{values}, fmt{fmt}
{}

const IO::Plugin* IO::Parameter::getFmt() const
{
    return fmt;
}

void IO::to_json(nlohmann::json& j, const Parameter& p)
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

void IO::from_json(const nlohmann::json& j, Parameter& p)
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
