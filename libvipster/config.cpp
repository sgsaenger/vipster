#include <config.h>
#include <json.hpp>
#include <fstream>
#include <locale>

using json = nlohmann::json;
using namespace Vipster;

PseMap Vipster::readPse()
{
    std::ifstream user_file{user_config};
    PseMap temp{true};
    if(user_file){
        json loc_file;
        loc_file << user_file;
        for(auto it=loc_file["PSE"].begin();it!=loc_file["PSE"].end();++it)
        {
            auto v = it.value();
            temp.emplace(it.key(),PseEntry{v["PWPP"],
                    v["CPPP"],v["CPNL"],v["Z"],
                    v["m"],v["bondcut"],v["covr"],
                    v["vdwr"],v["col"]});
        }
    }
    if(temp.find("X")==temp.end()){
        temp.emplace("X", PseEntry{"","","",0,0.,1.46,1.46,3.21,{{0.,0.,0.,1.}}});
    }
    return temp;
}

PseEntry& PseMap::operator [](const std::string& k)
{
    try{
        return at(k);
    }catch(const std::out_of_range& ex){
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
        emplace(k,pse.at("X"));
        return at(k);
    }
}
