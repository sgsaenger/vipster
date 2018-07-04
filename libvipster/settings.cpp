#include "settings.h"

using namespace Vipster;

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
