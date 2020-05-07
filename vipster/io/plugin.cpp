#include "plugin.h"
#include "plugins/xyz.h"
#include "plugins/pwinput.h"
#include "plugins/pwoutput.h"
#include "plugins/lmpinput.h"
#include "plugins/lmptrajec.h"
#include "plugins/cpmdinput.h"
#include "plugins/cube.h"
#include "plugins/xsf.h"
#include "plugins/orca.h"
#include "plugins/poscar.h"

using namespace Vipster;

IO::Plugins IO::defaultPlugins()
{
    return {
            &IO::XYZ,
            &IO::PWInput,
            &IO::PWOutput,
            &IO::LmpInput,
            &IO::LmpTrajec,
            &IO::CPInput,
            &IO::Cube,
            &IO::XSF,
            &IO::OrcaInput,
            &IO::Poscar
    };
}
