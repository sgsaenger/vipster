#include "preset.h"

using namespace Vipster;

IO::LmpPreset::LmpPreset(AtomStyle style,
                         bool bonds, bool angles, bool dihedrals, bool impropers)
    : style{style},
      bonds{bonds}, angles{angles},
      dihedrals{dihedrals}, impropers{impropers}
{}

IOFmt IO::LmpPreset::getFmt() const
{
    return IOFmt::LMP;
}

std::unique_ptr<IO::BasePreset> IO::LmpPreset::copy() const
{
    return std::make_unique<IO::LmpPreset>(*this);
}

void IO::LmpPreset::parseJson(const nlohmann::json& j)
{
    angles = j.value("angles", false);
    bonds = j.value("bonds", false);
    dihedrals = j.value("dihedrals", false);
    impropers = j.value("impropers", false);
    auto tmp = j.find("style");
    if(tmp != j.end()){
        auto stylestr = tmp->get<std::string>();
        for(const auto& pair: IO::LmpPreset::fmt2str){
            if(pair.second == stylestr){
                style = pair.first;
            }
        }
    }else{
        style = IO::LmpPreset::AtomStyle::Atomic;
    }
}

nlohmann::json IO::LmpPreset::toJson() const
{
    nlohmann::json j;
    j["angles"] = angles;
    j["bonds"] = bonds;
    j["dihedrals"] = dihedrals;
    j["impropers"] = impropers;
    j["style"] = IO::LmpPreset::fmt2str.at(style);
    return j;
}

const std::map<IO::LmpPreset::AtomStyle, std::string> IO::LmpPreset::fmt2str{
    {AtomStyle::Angle, "angle"},
    {AtomStyle::Atomic, "atomic"},
    {AtomStyle::Bond, "bond"},
    {AtomStyle::Charge, "charge"},
    {AtomStyle::Full, "full"},
    {AtomStyle::Molecular, "molecular"},
};
