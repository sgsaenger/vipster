#include "plugin.h"

using namespace Vipster;

IO::XYZConfig::XYZConfig(std::string n, Mode m, Data d)
    : BaseConfig{n}, filemode{m}, atomdata{d}
{}

IOFmt IO::XYZConfig::getFmt() const
{
    return IOFmt::XYZ;
}

std::unique_ptr<IO::BaseConfig> IO::XYZConfig::copy() const
{
    return std::make_unique<IO::XYZConfig>(*this);
}

void Vipster::IO::to_json(nlohmann::json& j, const IO::XYZConfig& c)
{
    j["filemode"] = c.filemode;
    j["atomdata"] = c.atomdata;
}

void Vipster::IO::from_json(const nlohmann::json& j, IO::XYZConfig& c)
{
    c.filemode = j.value("filemode", IO::XYZConfig::Mode::Step);
    c.atomdata = j.value("atomdata", IO::XYZConfig::Data::None);
}

void IO::XYZConfig::parseJson(const nlohmann::json::iterator& it)
{
    name = it.key();
    from_json(it.value(), *this);
}

nlohmann::json IO::XYZConfig::toJson() const
{
    nlohmann::json j;
    to_json(j, *this);
    return j;
}
