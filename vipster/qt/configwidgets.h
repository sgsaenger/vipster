#ifndef CONFIGWIDGETS_H
#define CONFIGWIDGETS_H

#include "configwidgets/pwconfig.h"
#include "configwidgets/lmpconfig.h"
#include "configwidgets/xyzconfig.h"
#include "configwidgets/cpconfig.h"
#include "configwidgets/poscarconfig.h"

inline std::map<Vipster::IOFmt, ConfigBase*> makeConfigWidgets()
{
    return {
        {Vipster::IOFmt::PWI, new PWConfig()},
        {Vipster::IOFmt::LMP, new LmpConfig()},
        {Vipster::IOFmt::XYZ, new XYZConfig()},
        {Vipster::IOFmt::CPI, new CPConfig()},
        {Vipster::IOFmt::POSCAR, new PoscarConfig()},
    };
}

#endif // CONFIGWIDGETS_H
