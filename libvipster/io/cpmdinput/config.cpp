#include "config.h"

using namespace Vipster;

IO::CPConfig::CPConfig(std::string name, bool angstrom, Scale scale)
    : BaseConfig{name}, angstrom{angstrom}, scale{scale}
{}

std::unique_ptr<IO::BaseConfig> IO::CPConfig::copy()
{
    return std::make_unique<IO::CPConfig>(*this);
}

void IO::CPConfig::parseJson(const nlohmann::json& j)
{
    from_json(j, *this);
}

nlohmann::json IO::CPConfig::toJson()
{
    nlohmann::json j;
    to_json(j, *this);
    return j;
}

void IO::from_json(const nlohmann::json& j, IO::CPConfig& c)
{
    c.name = j.at("name");
    c.angstrom = j.at("angstrom");
    c.scale = j.at("scale");
}

void IO::to_json(nlohmann::json& j, const IO::CPConfig& c)
{
    j["name"] = c.name;
    j["angstrom"] = c.angstrom;
    j["scale"] = c.scale;
}
