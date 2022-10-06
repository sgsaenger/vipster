#ifndef Settings_H
#define Settings_H

#include "global.h"

#include <string>

namespace Vipster{

template<typename T>
struct Setting{
    using ValueType = T;
    std::string name;
    T val;
};

struct Settings{
    Setting<bool>   overlap{"Check for overlapping atoms", true};
    Setting<bool>   atRadVdW{"Atom radius VdW", false};
    Setting<double> atRadFac{"Atom radius factor", bohrrad};
    Setting<double> bondRad{"Bond radius", bohrrad};
    Setting<bool>   showCell{"Show cell", true};
    Setting<bool>   antialias{"Antialiasing", true};
    Setting<bool>   perspective{"Perspective projection", false};
    Setting<bool>   rotCom{"Rotate around center of mass", false};
    Setting<size_t> animstep{"Animation step (ms)", 100};
    Setting<ColVec> bgCol{"Background color", ColVec{255, 255, 255, 255}};
    Setting<ColVec> selCol{"Selection color", ColVec{0, 0, 80, 80}};
    Setting<ColVec> milCol{"Miller-plane color", ColVec{130, 0, 0, 80}};
    Setting<ColVec> posCol{"Positive-isovalue color", ColVec{255, 0, 0, 155}};
    Setting<ColVec> negCol{"Negative-isovalue color", ColVec{0, 0, 255, 155}};
};

inline const Settings settings{};

}

#endif // Settings_H
