#include "config.h"

using namespace Vipster;

IO::CPConfig::CPConfig(std::string name, CPConfig::AtomFmt fmt)
    : BaseConfig{name}, fmt{fmt}
{}

IOFmt IO::CPConfig::getFmt() const
{
    return IOFmt::CPI;
}

std::unique_ptr<IO::BaseConfig> IO::CPConfig::copy() const
{
    return std::make_unique<IO::CPConfig>(*this);
}

void IO::CPConfig::parseJson(const nlohmann::json::iterator& it)
{
    name = it.key();
    from_json(it.value(), *this);
}

nlohmann::json IO::CPConfig::toJson() const
{
    nlohmann::json j;
    to_json(j, *this);
    return j;
}

void IO::from_json(const nlohmann::json& j, IO::CPConfig& c)
{
    c.fmt = j.value("fmt", IO::CPConfig::AtomFmt::Current);
}

void IO::to_json(nlohmann::json& j, const IO::CPConfig& c)
{
    j["fmt"] = c.fmt;
}
