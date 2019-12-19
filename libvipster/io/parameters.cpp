#include "parameters.h"
#include "../fileio.h"

using namespace Vipster;

IO::BaseParam::BaseParam(const struct Plugin* fmt,
                         CustomMap<std::string, ParamValue> &&values)
    : CustomMap{values}, fmt{fmt}
{}

const IO::Plugin* IO::BaseParam::getFmt() const
{
    return fmt;
}

void IO::to_json(nlohmann::json& j, const BaseParam& p)
{
    for(const auto &v: p){
        switch (v.second.index()) {
        default:
            break;
        }
    }
}

void IO::from_json(const nlohmann::json& j, BaseParam& p)
{
    for(auto &v: p){
        switch (v.second.index()) {
        default:
            break;
        }
    }
}
