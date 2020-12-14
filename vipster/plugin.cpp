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
#include "plugins/json.h"

namespace Vipster{

PluginList defaultPlugins()
{
    return {
            &Plugins::XYZ,
            &Plugins::PWInput,
            &Plugins::PWOutput,
            &Plugins::LmpInput,
            &Plugins::LmpTrajec,
            &Plugins::CPInput,
            &Plugins::Cube,
            &Plugins::XSF,
            &Plugins::OrcaInput,
            &Plugins::Poscar,
            &Plugins::JSON
    };
}
}
