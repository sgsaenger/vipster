#ifndef CONFIGWIDGETS_H
#define CONFIGWIDGETS_H

#include "configwidgets/pwconfig.h"
#include "configwidgets/lmpconfig.h"
#include "configwidgets/xyzconfig.h"
#include "configwidgets/cpconfig.h"
#include "configwidgets/poscarconfig.h"

inline std::map<const Vipster::IO::Plugin*, ConfigBase*> makeConfigWidgets()
{
    return {
        {&Vipster::IO::PWInput, new PWConfig()},
        {&Vipster::IO::LmpInput, new LmpConfig()},
        {&Vipster::IO::XYZ, new XYZConfig()},
        {&Vipster::IO::CPInput, new CPConfig()},
        {&Vipster::IO::Poscar, new PoscarConfig()},
    };
}

#endif // CONFIGWIDGETS_H
