#ifndef PSE_H
#define PSE_H

#include "global.h"
#include "json.hpp"

#include <string>
#include <map>

namespace Vipster{

struct Element{
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

void to_json(nlohmann::json& j,const Element& p);
void from_json(const nlohmann::json& j, Element& p);

class PeriodicTable:private std::map<std::string,Element>
{
    using map_t = std::map<std::string, Element>;
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
    PeriodicTable(std::initializer_list<PeriodicTable::value_type> il={}, bool r=false);
    Element& operator [](const std::string &k);
    iterator find_or_fallback(const std::string &k);
private:
    bool root;
};

void to_json(nlohmann::json& j,const PeriodicTable& p);
void from_json(const nlohmann::json& j, PeriodicTable& p);

extern PeriodicTable pte;

}

#endif // PSE_H
