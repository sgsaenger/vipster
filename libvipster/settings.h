#ifndef Settings_H
#define Settings_H

#include "json.hpp"
#include "global.h"
#include "bond.h"

#include <string>

namespace Vipster{

template<typename T>
struct Setting{
    using ValueType = T;
    std::string name;
    T val;
};

struct Settings{
    Setting<bool> atRadVdW{"Atom radius VdW", false};
    Setting<float> atRadFac{"Atom radius factor", bohrrad};
    Setting<float> bondRad{"Bond radius", bohrrad};
    Setting<float> bondCutFac{"Bond cutoff factor", 1.1f};
    Setting<BondFrequency> bondFreq{"Bond frequency", BondFrequency::Always};
    Setting<BondLevel> bondLvl{"Bond level", BondLevel::Cell};
    Setting<bool> showBonds{"Show bonds", true};
    Setting<bool> showCell{"Show cell", true};
    Setting<bool> antialias{"Antialiasing", true};
    Setting<bool> perspective{"Perspective projection", false};
    Setting<bool> rotCom{"Rotate around center of mass", false};
    Setting<size_t> animstep{"Animation step (ms)", 100};
    Setting<ColVec> selCol{"Selection color", ColVec{0, 0, 80, 80}};
    Setting<ColVec> milCol{"Miller-plane color", ColVec{130, 0, 0, 80}};
    Setting<ColVec> posCol{"Positive-isovalue color", ColVec{255, 0, 0, 155}};
    Setting<ColVec> negCol{"Negative-isovalue color", ColVec{0, 0, 255, 155}};
    Setting<std::string> PWPP{"Default PWScf PP-suffix", ""};
    Setting<std::string> CPPP{"Default CPMD PP-suffix", ""};
    Setting<std::string> CPNL{"Default CPMD Nonlocality", "LMAX=F"};
};

void to_json(nlohmann::json& j, const Settings& s);
void from_json(const nlohmann::json& j, Settings& s);

extern Settings settings;

}

#endif // Settings_H
