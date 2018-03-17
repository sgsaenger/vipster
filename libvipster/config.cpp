#include "config.h"
#include "json.hpp"
#include <fstream>
#include <locale>

using json = nlohmann::json;
using namespace Vipster;

namespace Vipster{
    extern const PseMap pse = readPse();
}

PseMap Vipster::readPse()
{
    PseMap temp{true};
    std::ifstream conf_file{user_config};
    if(!conf_file){
        conf_file = std::ifstream{sys_config};
    }
    if(conf_file){
        json loc_file;
        conf_file >> loc_file;
        for(auto it=loc_file["PSE"].begin();it!=loc_file["PSE"].end();++it)
        {
            auto v = it.value();
            temp.emplace(it.key(),PseEntry{v["PWPP"],
                    v["CPPP"],v["CPNL"],v["Z"],
                    v["m"],v["bondcut"],v["covr"],
                    v["vdwr"],v["col"]});
        }
    }
    // ensure fallback-value is present
    if(temp.find("")==temp.end()){
        temp.emplace("", PseEntry{"","","",0,0,0,1.46f,3.21f,{{0,0,0,255}}});
    }
    return temp;
}

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
                if(std::islower(test[0])){
                    test[0] = std::toupper(test[0],loc);
                }
                if(pse.find(test)!=pse.end())
                {
                    emplace(k,pse.at(test));
                    return at(k);
                }
            }
        }
        emplace(k,pse.at(""));
        return at(k);
    }
}
