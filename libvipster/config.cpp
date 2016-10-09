#include <config.h>
#include <json.hpp>
#include <fstream>
#include <locale>

using json = nlohmann::json;
using namespace Vipster;

std::map<std::string, PseEntry> Vipster::readPse()
{
    std::ifstream user_file{user_config};
    std::map<std::string,PseEntry> temp;
    if(user_file){
        json loc;
        loc << user_file;
        for(auto it=loc["PSE"].begin();it!=loc["PSE"].end();++it)
        {
            auto v = it.value();
            temp[it.key()]=PseEntry{v["PWPP"],
                    v["CPPP"],v["CPNL"],v["Z"],
                    v["m"],v["bondcut"],v["covr"],
                    v["vdwr"],v["col"]};
        }
    }else{
        temp["X"]=PseEntry{"","","",0,0.,1.46,1.46,3.21,{{0.,0.,0.,}}};
    }
    return temp;
}

PseEntry& PseMap::operator [](const std::string& k)
{
    try{
        return at(k);
    }catch(const std::out_of_range& ex){
        std::locale loc;
        for(size_t i=k.length();i>0;--i)
        {
            std::string test = k.substr(0,i);
            if(std::islower(test[0])){
                test[0] = std::toupper(test[0],loc);
            }
            if(internal->find(test)!=internal->end())
            {
                emplace(k,internal->at(test));
                return at(k);
            }
        }
        emplace(k,internal->at("X"));
        return at(k);
    }
}
