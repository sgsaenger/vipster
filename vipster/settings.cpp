#include "settings.h"
#include <nlohmann/json.hpp>

using namespace Vipster;

template<typename T>
void readSetting(const nlohmann::json& j, Setting<T>& s)
{
    auto tmp = j.find(s.name);
    if(tmp != j.end()){
        s.val = *tmp;
    }
}

namespace Vipster{
void from_json(const nlohmann::json& j, Settings& s){
    readSetting(j, s.overlap);
    readSetting(j, s.atRadFac);
    readSetting(j, s.atRadVdW);
    readSetting(j, s.bondRad);
    readSetting(j, s.showCell);
    readSetting(j, s.antialias);
    readSetting(j, s.perspective);
    readSetting(j, s.rotCom);
    readSetting(j, s.animstep);
    readSetting(j, s.bgCol);
    readSetting(j, s.selCol);
    readSetting(j, s.milCol);
    readSetting(j, s.posCol);
    readSetting(j, s.negCol);
}

void to_json(nlohmann::json& j, const Settings& s){
    j[s.overlap.name] = s.overlap.val;
    j[s.atRadFac.name] = s.atRadFac.val;
    j[s.atRadVdW.name] = s.atRadVdW.val;
    j[s.bondRad.name] = s.bondRad.val;
    j[s.showCell.name] = s.showCell.val;
    j[s.antialias.name] = s.antialias.val;
    j[s.perspective.name] = s.perspective.val;
    j[s.rotCom.name] = s.rotCom.val;
    j[s.animstep.name] = s.animstep.val;
    j[s.bgCol.name] = s.bgCol.val;
    j[s.selCol.name] = s.selCol.val;
    j[s.milCol.name] = s.milCol.val;
    j[s.posCol.name] = s.posCol.val;
    j[s.negCol.name] = s.negCol.val;
}
}
