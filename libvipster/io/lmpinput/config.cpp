#include "config.h"

using namespace Vipster;

IO::LmpConfig::LmpConfig(std::string name, AtomStyle style,
                         bool bonds, bool angles, bool dihedrals, bool impropers)
    : BaseConfig{name},
      style{style},
      bonds{bonds}, angles{angles},
      dihedrals{dihedrals}, impropers{impropers}
{}

IOFmt IO::LmpConfig::getFmt() const
{
    return IOFmt::LMP;
}

std::unique_ptr<IO::BaseConfig> IO::LmpConfig::copy() const
{
    return std::make_unique<IO::LmpConfig>(*this);
}

void IO::LmpConfig::parseJson(const nlohmann::json::iterator& it)
{
    name = it.key();
    from_json(it.value(), *this);
}

nlohmann::json IO::LmpConfig::toJson() const
{
    nlohmann::json j;
    to_json(j, *this);
    return j;
}

void IO::from_json(const nlohmann::json& j, IO::LmpConfig& c)
{
    c.angles = j.value("angles", false);
    c.bonds = j.value("bonds", false);
    c.dihedrals = j.value("dihedrals", false);
    c.impropers = j.value("impropers", false);
    auto style = j.find("style");
    if(style != j.end()){
        auto stylestr = style->get<std::string>();
        for(const auto& pair: IO::LmpConfig::fmt2str){
            if(pair.second == stylestr){
                c.style = pair.first;
            }
        }
    }else{
        c.style = IO::LmpConfig::AtomStyle::Atomic;
    }
}

void IO::to_json(nlohmann::json &j, const IO::LmpConfig& c)
{
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
