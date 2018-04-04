#include "config.h"
#include "iowrapper.h"
#include <fstream>
#include <locale>
#include <tuple>

using json = nlohmann::json;
using namespace Vipster;

bool readConfig();

namespace Vipster{
PseMap pse;
//Settings settings;
Parameters params;
bool config_loaded = readConfig();
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
                if(std::islower(test[0]) != 0){
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

bool readConfig()
{
    PseMap pse{true};
//    Settings settings;
    Parameters params;
    std::ifstream conf_file{user_config};
    if(!conf_file){
        conf_file = std::ifstream{sys_config};
    }
    if(conf_file){
        json loc_file;
        conf_file >> loc_file;
        // PSE
        for(auto it=loc_file["PSE"].begin();it!=loc_file["PSE"].end();++it)
        {
            auto v = it.value();
            pse.emplace(it.key(),PseEntry{v["PWPP"],
                    v["CPPP"],v["CPNL"],v["Z"],
                    v["m"],v["bondcut"],v["covr"],
                    v["vdwr"],v["col"]});
        }
        // General settings
        // TODO: figure out c++17-support on CI
//        for(auto it=loc_file["Settings"].begin();it!=loc_file["Settings"].end();++it)
//        {
//            if(it.value().is_boolean()){
//                settings.insert({it.key(), it.value().get<bool>()});
//            }else if(it.value().is_number_float()){
//                settings.insert({it.key(), it.value().get<float>()});
//            }else if(it.value().is_number_integer()){
//                settings.insert({it.key(), it.value().get<int>()});
//            }else{
//                settings.insert({it.key(), it.value().get<std::string>()});
//            }
//        }
        // Parameter sets
        for(const auto& pw:loc_file.at("Parameters").at("PWI")){
            params.emplace(std::make_pair(IOFmt::PWI, std::make_unique<IO::PWParam>(pw.get<IO::PWParam>())));
        }
    }
    // ensure fallback-value is present
    if(pse.find("")==pse.end()){
        pse.emplace("", PseEntry{"","","",0,0,0,1.46f,3.21f,{{0,0,0,255}}});
    }
    Vipster::pse = std::move(pse);
//    Vipster::settings = std::move(settings);
    Vipster::params = std::move(params);
    return true;
}
