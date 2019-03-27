#include "config.h"

using namespace Vipster;

IO::PWConfig::PWConfig(std::string name, PWConfig::AtomFmt atoms, CellFmt cell)
    : BaseConfig{name}, atoms{atoms}, cell{cell}
{}

IOFmt IO::PWConfig::getFmt() const
{
    return IOFmt::PWI;
}

std::unique_ptr<IO::BaseConfig> IO::PWConfig::copy() const
{
    return std::make_unique<IO::PWConfig>(*this);
}

void Vipster::IO::to_json(nlohmann::json &j, const IO::PWConfig& c)
{
    j["atomfmt"] = c.atoms;
    j["cellfmt"] = c.cell;
}

void Vipster::IO::from_json(const nlohmann::json &j, IO::PWConfig& c)
{
    c.atoms = j.value("atomfmt", IO::PWConfig::AtomFmt::Current);
    c.cell = j.value("cellfmt", IO::PWConfig::CellFmt::Current);
}

void IO::PWConfig::parseJson(const nlohmann::json::iterator& it)
{
    name = it.key();
    from_json(it.value(), *this);
}

nlohmann::json IO::PWConfig::toJson() const
{
    nlohmann::json j;
    to_json(j, *this);
    return j;
}
