#include "preset.h"
#include "plugin.h"

using namespace Vipster;

IO::PWPreset::PWPreset(PWPreset::AtomFmt atoms, CellFmt cell)
    : atoms{atoms}, cell{cell}
{}

const IO::Plugin* IO::PWPreset::getFmt() const
{
    return &PWInput;
}

std::unique_ptr<IO::BasePreset> IO::PWPreset::copy() const
{
    return std::make_unique<IO::PWPreset>(*this);
}

void IO::PWPreset::parseJson(const nlohmann::json& j)
{
    atoms = j.value("atomfmt", IO::PWPreset::AtomFmt::Active);
    cell = j.value("cellfmt", IO::PWPreset::CellFmt::Active);
}

nlohmann::json IO::PWPreset::toJson() const
{
    nlohmann::json j;
    j["atomfmt"] = atoms;
    j["cellfmt"] = cell;
    return j;
}
