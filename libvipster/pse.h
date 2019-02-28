#ifndef PSE_H
#define PSE_H

#include "global.h"
#include "json.hpp"

#include <string>
#include <map>

namespace Vipster{

struct PseEntry{
    std::string     PWPP;
    std::string     CPPP;
    std::string     CPNL;
    unsigned int    Z;
    float           m;
    float           bondcut;
    float           covr;
    float           vdwr;
    ColVec          col;
};

void to_json(nlohmann::json& j,const PseEntry& p);
void from_json(const nlohmann::json& j, PseEntry& p);

class PseMap:private std::map<std::string,PseEntry>
{
    using map_t = std::map<std::string, PseEntry>;
public:
    using map_t::begin;
    using map_t::end;
    using map_t::at;
    using map_t::emplace;
    using map_t::find;
    using map_t::size;
    using map_t::key_type;
    using map_t::mapped_type;
    using map_t::value_type;
    using map_t::insert_or_assign;
    PseMap(std::initializer_list<PseMap::value_type> il={}, bool r=false);
    PseEntry& operator [](const std::string &k);
private:
    bool root;
};

void to_json(nlohmann::json& j,const PseMap& p);
void from_json(const nlohmann::json& j, PseMap& p);

extern PseMap pse;

}

#endif // PSE_H
