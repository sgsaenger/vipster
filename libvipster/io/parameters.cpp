#include "parameters.h"
#include "io/plugin.h"
#include "../fileio.h"

using namespace Vipster;

void IO::to_json(nlohmann::json& j, const BaseParam& p)
{
    j = p.toJson();
}

void IO::from_json(const nlohmann::json& j, BaseParam& p)
{
    p.parseJson(j);
}
