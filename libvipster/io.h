#ifndef IOWRAPPER
#define IOWRAPPER

#include <map>

#include "molecule.h"
#include "io/fmt.h"
#include "io/plugin.h"
#include "io/xyz/plugin.h"
#include "io/pwinput/plugin.h"
#include "io/pwoutput/plugin.h"
#include "io/lmpinput/plugin.h"
#include "io/lmptrajec/plugin.h"
#include "io/cpmdinput/plugin.h"
#include "io/cube/plugin.h"
#include "io/xsf/plugin.h"
#include "io/orca/plugin.h"

//TODO: check std::ios_base::sync_with_stdio(false)
namespace Vipster{
    const std::map<IOFmt, IO::Plugin const *const> IOPlugins{
            {IOFmt::XYZ, &IO::XYZ},
            {IOFmt::PWI, &IO::PWInput},
            {IOFmt::PWO, &IO::PWOutput},
            {IOFmt::LMP, &IO::LmpInput},
            {IOFmt::DMP, &IO::LmpTrajec},
            {IOFmt::CPI, &IO::CPInput},
            {IOFmt::CUBE, &IO::Cube},
            {IOFmt::XSF, &IO::XSF},
            {IOFmt::ORCA, &IO::OrcaInput}
    };
    IO::Data readFile(const std::string &fn, IOFmt fmt);
    IO::Data readFile(const std::string &fn, IOFmt fmt, std::string name);
    bool     writeFile(const std::string &fn, IOFmt fmt, const Molecule &m,
                       const IO::BaseParam *p=nullptr,
                       const IO::BaseConfig *c=nullptr,
                       IO::State state={});
}

#endif // IOWRAPPER

