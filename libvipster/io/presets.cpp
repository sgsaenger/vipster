#include "presets.h"
#include "../fileio.h"

using namespace Vipster;

void IO::to_json(nlohmann::json& j, const BasePreset& p)
{
    j = p.toJson();
}

void IO::from_json(const nlohmann::json& j, BasePreset& p)
{
    p.parseJson(j);
}
