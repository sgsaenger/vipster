#ifndef CONFIGWIDGETS_H
#define CONFIGWIDGETS_H

#include "presetwidgets/pwpreset.h"
#include "presetwidgets/lmppreset.h"
#include "presetwidgets/xyzpreset.h"
#include "presetwidgets/cppreset.h"
#include "presetwidgets/poscarpreset.h"

inline std::map<const Vipster::IO::Plugin*, PresetBase*> makePresetWidgets()
{
    return {
        {&Vipster::IO::PWInput, new PWPreset()},
        {&Vipster::IO::LmpInput, new LmpPreset()},
        {&Vipster::IO::XYZ, new XYZPreset()},
        {&Vipster::IO::CPInput, new CPPreset()},
        {&Vipster::IO::Poscar, new PoscarPreset()},
    };
}

#endif // CONFIGWIDGETS_H
