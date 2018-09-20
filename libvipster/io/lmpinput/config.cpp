#include "config.h"

using namespace Vipster;

IO::LmpConfig::LmpConfig(std::string name, AtomStyle style,
                         bool bonds, bool angles, bool dihedrals, bool impropers)
    : BaseConfig{name},
      style{style},
      bonds{bonds}, angles{angles},
      dihedrals{dihedrals}, impropers{impropers}
{}

std::unique_ptr<IO::BaseConfig> IO::LmpConfig::copy()
{
    return std::make_unique<IO::LmpConfig>(*this);
}

void IO::LmpConfig::parseJson(const nlohmann::json& j)
{
    from_json(j, *this);
}

nlohmann::json IO::LmpConfig::toJson()
{
    nlohmann::json j;
    to_json(j, *this);
    return j;
}

void IO::from_json(const nlohmann::json& j, IO::LmpConfig& c)
{
    c.name = j.at("name");
    c.angles = j.at("angles");
    c.bonds = j.at("bonds");
    c.dihedrals = j.at("dihedrals");
    c.impropers = j.at("impropers");
    for(const auto& pair: IO::LmpConfig::fmt2str){
        if(pair.second == j.at("style")){
            c.style = pair.first;
        }
    }
}

void IO::to_json(nlohmann::json &j, const IO::LmpConfig& c)
{
    j["name"] = c.name;
    j["angles"] = c.angles;
    j["bonds"] = c.bonds;
    j["dihedrals"] = c.dihedrals;
    j["impropers"] = c.impropers;
    j["style"] = IO::LmpConfig::fmt2str.at(c.style);
}

const std::map<IO::LmpConfig::AtomStyle, std::string> IO::LmpConfig::fmt2str{
    {AtomStyle::Angle, "angle"},
    {AtomStyle::Atomic, "atomic"},
    {AtomStyle::Bond, "bond"},
    {AtomStyle::Charge, "charge"},
    {AtomStyle::Full, "full"},
    {AtomStyle::Molecular, "molecular"},
};
