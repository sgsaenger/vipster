#include "periodictable.h"
#include <cctype>
#include <optional>
#include <nlohmann/json.hpp>

using namespace Vipster;

bool Vipster::operator==(const Element &lhs, const Element &rhs)
{
    return std::tie(lhs.PWPP, lhs.CPPP, lhs.CPNL, lhs.Z, lhs.m,
                    lhs.bondcut, lhs.covr, lhs.vdwr, lhs.col)
            ==
           std::tie(rhs.PWPP, rhs.CPPP, rhs.CPNL, rhs.Z, rhs.m,
                    rhs.bondcut, rhs.covr, rhs.vdwr, rhs.col);
}

PeriodicTable::PeriodicTable(std::initializer_list<PeriodicTable::value_type> il,
                             const PeriodicTable *r)
    : StaticMap{il}, root{r}
{}

PeriodicTable::iterator PeriodicTable::find_or_fallback(const std::string &k)
{
    // send lookup either to specific root or hard-coded fallback-table
    const PeriodicTable& root = this->root? *this->root : Vipster::pte;
    auto entry = find(k);
    if(entry != end()){
        return entry;
    }else{
        // if key is ONLY a number, interpret as atomic number
        char *p;
        std::size_t Z = std::strtoul(k.c_str(), &p, 10);
        if(*p){
            // found a derived/custom name, try to guess base-name
            bool islower = std::islower(k[0]);
            // gradually ignore appended letters until we reach a matching atom type
            for(size_t i=k.length(); i>0; --i)
            {
                auto search = [&](const std::string& key)->std::optional<const Element>{
                    const_iterator tmp = find(key);
                    // local lookup
                    if(tmp != end()){
                        return {tmp->second};
                    }
                    // lookup in fallback table
                    tmp = root.find(key);
                    if(tmp != root.end()){
                        return {tmp->second};
                    }
                    return {};
                };
                if(islower){
                    // first check values starting with a uppercase letter to match defaults
                    std::string tmp = k.substr(0,i);
                    tmp[0] = std::toupper(k[0]);
                    auto test = search(tmp);
                    if(test){
                        return emplace(k, test.value()).first;
                    }
                }
                std::string tmp = k.substr(0,i);
                auto test = search(tmp);
                if(test){
                    return emplace(k, test.value()).first;
                }
            }
        }else{
            // interpret atomic number
            for(const auto& pair: root){
                if(pair.second.Z == Z){
                    return emplace(k, pair.second).first;
                }
            }
        }
    }
    return emplace(k, root.at("")).first;
}

Element& PeriodicTable::operator [](const std::string& k)
{
    return find_or_fallback(k)->second;
}

namespace Vipster{
void to_json(nlohmann::json& j, const Element& p)
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

void from_json(const nlohmann::json& j, Element& p)
{
    p.PWPP = j.value("PWPP", "");
    p.CPPP = j.value("CPPP", "");
    p.CPNL = j.value("CPNL", "");
    p.Z = j.value("Z", 0u);
    p.m = j.value("m", 0.0);
    p.bondcut = j.value("bondcut", -1.);
    p.covr = j.value("covr", 1.46);
    p.vdwr = j.value("vdwr", 3.21);
    p.col = j.value("col", ColVec{0,0,0,255});
}

void to_json(nlohmann::json& j, const PeriodicTable& m)
{
    for(const auto& e: m){
        j[e.first] = e.second;
    }
}

void from_json(const nlohmann::json& j, PeriodicTable& m)
{
    for(auto it=j.begin(); it!=j.end(); ++it){
        m.insert_or_assign(it.key(), it.value());
    }
}
}
