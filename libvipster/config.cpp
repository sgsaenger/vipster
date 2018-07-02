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
Settings settings;
Parameters params;
Configs configs;
static bool config_loaded = readConfig();
}

BaseParam::BaseParam(std::string name)
    :name{name}
{}

BaseConfig::BaseConfig(std::string name)
    :name{name}
{}

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
    Settings settings;
    Parameters params;
    Configs configs;
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
        settings = loc_file["Settings"];
        // Parameter sets
        for(const auto& pw:loc_file["Parameters"]["PWI"]){
            params.emplace(IOFmt::PWI, std::make_unique<IO::PWParam>(pw.get<IO::PWParam>()));
        }
    }
    // ensure fallback-values are present
    if(pse.find("")==pse.end()){
        pse.emplace("", PseEntry{"","","",0,0,0,1.46f,3.21f,{{0,0,0,255}}});
    }
    if(configs.find(IOFmt::XYZ) == configs.end()){
        configs.emplace(IOFmt::XYZ, std::make_unique<IO::XYZConfig>(IO::XYZConfigDefault));
    }
    if(params.find(IOFmt::PWI) == params.end()){
        params.emplace(IOFmt::PWI, std::make_unique<IO::PWParam>(IO::PWParamDefault));
    }
    if(configs.find(IOFmt::PWI) == configs.end()){
        configs.emplace(IOFmt::PWI, std::make_unique<IO::PWConfig>(IO::PWConfigDefault));
    }
    if(configs.find(IOFmt::LMP) == configs.end()){
        configs.emplace(IOFmt::LMP, std::make_unique<IO::LmpConfig>(IO::LmpConfigDefault));
    }
    Vipster::pse = std::move(pse);
    Vipster::settings = std::move(settings);
    Vipster::params = std::move(params);
    Vipster::configs = std::move(configs);
    return true;
}

template<typename T>
void jsonToSetting(const nlohmann::json& j, T& s){
    using ValueType = typename T::ValueType;
    auto it = j.find(s.name);
    if(it != j.end()){
        s.val = it.value().template get<ValueType>();
    }
}

void Vipster::from_json(const nlohmann::json& j, Settings& s){
    jsonToSetting(j, s.atRadFac);
    jsonToSetting(j, s.atRadVdW);
    jsonToSetting(j, s.bondRad);
    jsonToSetting(j, s.bondCutFac);
    jsonToSetting(j, s.bondFreq);
    jsonToSetting(j, s.bondLvl);
    jsonToSetting(j, s.showBonds);
    jsonToSetting(j, s.showCell);
    jsonToSetting(j, s.antialias);
    jsonToSetting(j, s.perspective);
    jsonToSetting(j, s.selCol);
    jsonToSetting(j, s.PWPP);
    jsonToSetting(j, s.CPPP);
    jsonToSetting(j, s.CPNL);
}
