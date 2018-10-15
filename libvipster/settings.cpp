#include "settings.h"

using namespace Vipster;

Vipster::Settings Vipster::settings;

template<typename T>
void readSetting(const nlohmann::json& j, Setting<T>& s)
{
    auto tmp = j.find(s.name);
    if(tmp != j.end()){
        s.val = *tmp;
    }
}

void Vipster::from_json(const nlohmann::json& j, Settings& s){
    readSetting(j, s.atRadFac);
    readSetting(j, s.atRadVdW);
    readSetting(j, s.bondRad);
    readSetting(j, s.bondCutFac);
    readSetting(j, s.bondFreq);
    readSetting(j, s.bondLvl);
    readSetting(j, s.showBonds);
    readSetting(j, s.showCell);
    readSetting(j, s.antialias);
    readSetting(j, s.perspective);
    readSetting(j, s.rotCom);
    readSetting(j, s.animstep);
    readSetting(j, s.selCol);
    readSetting(j, s.milCol);
    readSetting(j, s.posCol);
    readSetting(j, s.negCol);
    readSetting(j, s.PWPP);
    readSetting(j, s.CPPP);
    readSetting(j, s.CPNL);
}

void Vipster::to_json(nlohmann::json& j, const Settings& s){
    j[s.atRadFac.name] = s.atRadFac.val;
    j[s.atRadVdW.name] = s.atRadVdW.val;
    j[s.bondRad.name] = s.bondRad.val;
    j[s.bondCutFac.name] = s.bondCutFac.val;
    j[s.bondFreq.name] = s.bondFreq.val;
    j[s.bondLvl.name] = s.bondLvl.val;
    j[s.showBonds.name] = s.showBonds.val;
    j[s.showCell.name] = s.showCell.val;
    j[s.antialias.name] = s.antialias.val;
    j[s.perspective.name] = s.perspective.val;
    j[s.rotCom.name] = s.rotCom.val;
    j[s.animstep.name] = s.animstep.val;
    j[s.selCol.name] = s.selCol.val;
    j[s.milCol.name] = s.milCol.val;
    j[s.posCol.name] = s.posCol.val;
    j[s.negCol.name] = s.negCol.val;
    j[s.PWPP.name] = s.PWPP.val;
    j[s.CPPP.name] = s.CPPP.val;
    j[s.CPNL.name] = s.CPNL.val;
}
