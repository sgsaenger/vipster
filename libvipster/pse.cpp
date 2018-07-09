#include "pse.h"
#include "configfile.h"
#include <locale>

using namespace Vipster;

PseEntry& PseMap::operator [](const std::string& k)
{
    try{
        return at(k);
    }catch(const std::out_of_range&){
        if(!root){
            std::locale loc;
            for(size_t i=k.length();i>0;--i)
            {
                std::string test = k.substr(0,i);
                if(std::islower(test[0]) != 0){
                    test[0] = std::toupper(test[0],loc);
                }
                if(Vipster::pse.find(test) != Vipster::pse.end())
                {
                    emplace(k, Vipster::pse.at(test));
                    return at(k);
                }
            }
        }
        emplace(k, Vipster::pse.at(""));
        return at(k);
    }
}

void Vipster::to_json(nlohmann::json& j, const PseEntry& p)
{
    j["PWPP"] = p.PWPP;
    j["CPPP"] = p.CPPP;
    j["CPNL"] = p.CPNL;
    j["Z"] = p.Z;
    j["m"] = p.m;
    j["bondcut"] = p.bondcut;
    j["covr"] = p.covr;
    j["vdwr"] = p.vdwr;
    j["col"] = p.col;
}

void Vipster::from_json(const nlohmann::json& j, PseEntry& p)
{
    p.PWPP = j.at("PWPP");
    p.CPPP = j.at("CPPP");
    p.CPNL = j.at("CPNL");
    p.Z = j.at("Z");
    p.m = j.at("m");
    p.bondcut = j.at("bondcut");
    p.covr = j.at("covr");
    p.vdwr = j.at("vdwr");
    p.col = j.at("col");
}

void Vipster::to_json(nlohmann::json& j, const PseMap& m)
{
    for(const auto& e: m){
        j[e.first] = e.second;
    }
}

void Vipster::from_json(const nlohmann::json& j, PseMap& m)
{
    for(auto it=j.begin(); it!=j.end(); ++it){
        m.emplace(it.key(), it.value());
    }
}
