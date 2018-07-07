#include "config.h"

using namespace Vipster;

IO::PWConfig::PWConfig(std::string name, PWConfig::AtomFmt atoms, CellFmt cell)
    : BaseConfig{name}, atoms{atoms}, cell{cell}
{}

std::unique_ptr<IO::BaseConfig> IO::PWConfig::copy()
{
    return std::make_unique<IO::PWConfig>(*this);
}

void Vipster::IO::to_json(nlohmann::json &j, const IO::PWConfig& c)
{
    j["name"] = c.name;
    j["atomfmt"] = c.atoms;
    j["cellfmt"] = c.cell;
}

void Vipster::IO::from_json(const nlohmann::json &j, IO::PWConfig& c)
{
    c.name = j.at("name");
    c.atoms = j.at("atomfmt");
    c.cell = j.at("cellfmt");
}

void IO::PWConfig::parseJson(const nlohmann::json& j)
{
    from_json(j, *this);
}

nlohmann::json IO::PWConfig::toJson()
{
    nlohmann::json j;
    to_json(j, *this);
    return j;
}
