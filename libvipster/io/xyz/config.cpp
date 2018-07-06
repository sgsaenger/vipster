#include "plugin.h"

using namespace Vipster;

IO::XYZConfig::XYZConfig(std::string n, Mode m, Data d)
    : BaseConfig{n}, filemode{m}, atomdata{d}
{}

std::unique_ptr<BaseConfig> IO::XYZConfig::copy()
{
    return std::make_unique<IO::XYZConfig>(*this);
}

void Vipster::IO::to_json(nlohmann::json& j, const IO::XYZConfig& c)
{
    j["name"] = c.name;
    j["filemode"] = c.filemode;
    j["atomdata"] = c.atomdata;
}

void Vipster::IO::from_json(const nlohmann::json& j, IO::XYZConfig& c)
{
    c.name = j.at("name");
    c.filemode = j.at("filemode");
    c.atomdata = j.at("atomdata");
}

void IO::XYZConfig::parseJson(const nlohmann::json& j)
{
    from_json(j, *this);
}
