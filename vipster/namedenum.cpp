#include "namedenum.h"
#include "global.h"

#include <algorithm>

using namespace Vipster;

NamedEnum::NamedEnum(int value, const std::vector<std::string> &names)
    : val{value}
{
    if(value > names.size()){
        throw Error{"NamedEnum value out of range"};
    }
    std::vector<std::pair<int, std::string>> tmp;
    for(size_t i=0; i<names.size(); ++i){
        tmp.emplace_back(i, names[i]);
    }
    StaticMap::operator=(StaticMap(tmp.begin(), tmp.end()));
}

NamedEnum::operator int() const
{
    return val;
}

int NamedEnum::value() const
{
    return val;
}

NamedEnum::operator const std::string&() const
{
    return at(val);
}

const std::string& NamedEnum::name() const
{
    return at(val);
}

std::vector<std::string> NamedEnum::names() const
{
    std::vector<std::string> n{};
    n.reserve(size());
    for(auto [k, v]: *this){
        n.push_back(v);
    }
    return n;
}

NamedEnum& NamedEnum::operator=(int v)
{
    if(v > size()){
        throw Error{"NamedEnum value out of range"};
    }
    val = v;
    return *this;
}

NamedEnum& NamedEnum::operator=(const std::string& s)
{
    if(auto pos = find_if(begin(), end(),
          [&](auto &p){return p.second == s;}); pos != end()){
        val = pos->first;
        return *this;
    }else{
        throw Error{"NamedEnum name unknown"};
    }
}
