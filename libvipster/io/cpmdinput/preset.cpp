#include "preset.h"

using namespace Vipster;

IO::CPPreset::CPPreset(CPPreset::AtomFmt fmt)
    : fmt{fmt}
{}

IOFmt IO::CPPreset::getFmt() const
{
    return IOFmt::CPI;
}

std::unique_ptr<IO::BasePreset> IO::CPPreset::copy() const
{
    return std::make_unique<IO::CPPreset>(*this);
}

void IO::CPPreset::parseJson(const nlohmann::json& j)
{
    fmt = j.value("fmt", IO::CPPreset::AtomFmt::Current);
}

nlohmann::json IO::CPPreset::toJson() const
{
    nlohmann::json j;
    j["fmt"] = fmt;
    return j;
}
