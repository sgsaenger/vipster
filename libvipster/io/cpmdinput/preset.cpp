#include "preset.h"
#include "plugin.h"

using namespace Vipster;

IO::CPPreset::CPPreset(CPPreset::AtomFmt fmt)
    : fmt{fmt}
{}

const IO::Plugin* IO::CPPreset::getFmt() const
{
    return &CPInput;
}

std::unique_ptr<IO::BasePreset> IO::CPPreset::copy() const
{
    return std::make_unique<IO::CPPreset>(*this);
}

void IO::CPPreset::parseJson(const nlohmann::json& j)
{
    fmt = j.value("fmt", IO::CPPreset::AtomFmt::Active);
}

nlohmann::json IO::CPPreset::toJson() const
{
    nlohmann::json j;
    j["fmt"] = fmt;
    return j;
}
