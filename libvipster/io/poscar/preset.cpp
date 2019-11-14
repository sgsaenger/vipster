#include "preset.h"
#include "plugin.h"

using namespace Vipster;

IO::PoscarPreset::PoscarPreset(bool s, bool c)
    : selective{s}, cartesian{c}
{}

const IO::Plugin *IO::PoscarPreset::getFmt() const
{
    return &Poscar;
}

std::unique_ptr<IO::BasePreset> IO::PoscarPreset::copy() const
{
    return std::make_unique<IO::PoscarPreset>(*this);
}

void IO::PoscarPreset::parseJson(const nlohmann::json &j)
{
    selective = j.value("selective", true);
    cartesian = j.value("cartesian", false);
}

nlohmann::json IO::PoscarPreset::toJson() const
{
    nlohmann::json j;
    j["selective"] = selective;
    j["cartesian"] = cartesian;
    return j;
}
