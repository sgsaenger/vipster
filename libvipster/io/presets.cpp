#include "presets.h"
#include "../fileio.h"

using namespace Vipster;

IO::CustomEnum::CustomEnum(int value, const std::vector<std::string> &names)
    : value{value}
{
    if(value > names.size()){
        throw Error{"CustomEnum value out of range"};
    }
    std::vector<std::pair<int, std::string>> tmp;
    for(size_t i=0; i<names.size(); ++i){
        tmp.emplace_back(i, names[i]);
    }
    CustomMap::operator=(CustomMap(tmp.begin(), tmp.end()));
}

IO::CustomEnum::operator int() const
{
    return value;
}

IO::CustomEnum::operator const std::string&() const
{
    return at(value);
}

IO::CustomEnum& IO::CustomEnum::operator=(int v)
{
    if(v > size()){
        throw IO::Error{"CustomEnum value out of range"};
    }
    value = v;
    return *this;
}

IO::CustomEnum& IO::CustomEnum::operator=(const std::string& s)
{
    if(auto pos = find_if(begin(), end(), [&](auto &p){return p.second == s;}); pos != end()){
        value = pos->first;
        return *this;
    }else{
        throw Error{"CustomEnum name unknown"};
    }
}

IO::Preset::Preset(const struct Plugin* fmt,
                           CustomMap<std::string, PresetValue> &&values)
    : CustomMap{values}, fmt{fmt}
{}

const IO::Plugin* IO::Preset::getFmt() const
{
    return fmt;
}

void IO::to_json(nlohmann::json& j, const Preset& p)
{
    for(const auto &v: p){
        switch (v.second.index()) {
        case Preset::i_bool:
            j[v.first] = std::get<Preset::i_bool>(v.second);
            break;
        case Preset::i_uint:
            j[v.first] = std::get<Preset::i_uint>(v.second);
            break;
        case Preset::i_enum:
            j[v.first] = static_cast<const std::string&>(std::get<Preset::i_enum>(v.second));
            break;
        default:
            break;
        }
    }
}

void IO::from_json(const nlohmann::json& j, Preset& p)
{
    for(auto &v: p){
        switch(v.second.index()) {
        case Preset::i_bool:
            v.second = j.value(v.first, std::get<Preset::i_bool>(v.second));
            break;
        case Preset::i_uint:
            v.second = j.value(v.first, std::get<Preset::i_uint>(v.second));
            break;
        case Preset::i_enum:
        {
            auto &val = std::get<Preset::i_enum>(v.second);
            if(auto pos = j.find(v.first); pos != j.end()){
                try{
                    val = pos->get<std::string>();
                }catch(IO::Error){
                    // ignore
                }
            }
        }
            break;
        default:
            break;
        }
    }
}
