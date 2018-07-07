#include "settings.h"

using namespace Vipster;

void Vipster::from_json(const nlohmann::json& j, Settings& s){
    s.atRadFac.val = j.at(s.atRadFac.name);
    s.atRadVdW.val = j.at(s.atRadVdW.name);
    s.bondRad.val = j.at(s.bondRad.name);
    s.bondCutFac.val = j.at(s.bondCutFac.name);
    s.bondFreq.val = j.at(s.bondFreq.name);
    s.bondLvl.val = j.at(s.bondLvl.name);
    s.showBonds.val = j.at(s.showBonds.name);
    s.showCell.val = j.at(s.showCell.name);
    s.antialias.val = j.at(s.antialias.name);
    s.perspective.val = j.at(s.perspective.name);
    s.selCol.val = j.at(s.selCol.name);
    s.PWPP.val = j.at(s.PWPP.name);
    s.CPPP.val = j.at(s.CPPP.name);
    s.CPNL.val = j.at(s.CPNL.name);
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
    j[s.selCol.name] = s.selCol.val;
    j[s.PWPP.name] = s.PWPP.val;
    j[s.CPPP.name] = s.CPPP.val;
    j[s.CPNL.name] = s.CPNL.val;
}
