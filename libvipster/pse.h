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
public:
    using std::map<std::string,PseEntry>::begin;
    using std::map<std::string,PseEntry>::end;
    using std::map<std::string,PseEntry>::at;
    using std::map<std::string,PseEntry>::emplace;
    using std::map<std::string,PseEntry>::find;
    using std::map<std::string,PseEntry>::size;
    PseMap(bool r=false):root(r){}
    PseEntry& operator [](const std::string &k);
private:
    bool root;
};

void to_json(nlohmann::json& j,const PseMap& p);
void from_json(const nlohmann::json& j, PseMap& p);

extern PseMap pse;

}

#endif // PSE_H
