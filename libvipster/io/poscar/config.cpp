#include "plugin.h"
#include "plugin.h"

using namespace Vipster;

IO::PoscarConfig::PoscarConfig(std::string n, bool s, bool c)
    : BaseConfig{n}, selective{s}, cartesian{c}
{}

const IO::Plugin *IO::PoscarConfig::getFmt() const
{
    return &Poscar;
}

std::unique_ptr<IO::BaseConfig> IO::PoscarConfig::copy() const
{
    return std::make_unique<IO::PoscarConfig>(*this);
}

void Vipster::IO::to_json(nlohmann::json &j, const IO::PoscarConfig &c)
{
    j["selective"] = c.selective;
    j["cartesian"] = c.cartesian;
}

void Vipster::IO::from_json(const nlohmann::json &j, IO::PoscarConfig &c)
{
    c.selective = j.value("selective", true);
    c.cartesian = j.value("cartesian", false);
}

void IO::PoscarConfig::parseJson(const nlohmann::json::iterator &it)
{
    name = it.key();
    from_json(it.value(), *this);
}

nlohmann::json IO::PoscarConfig::toJson() const
{
    nlohmann::json j;
    to_json(j, *this);
    return j;
}
