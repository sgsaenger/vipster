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
