#ifndef IOWRAPPER
#define IOWRAPPER

#include <map>

#include "molecule.h"
#include "io/fmt.h"
#include "io/plugin.h"
#include "io/xyz/xyz.h"
#include "io/pwinput/pwinput.h"
#include "io/pwoutput/pwoutput.h"
#include "io/lmpinput/lmpinput.h"
#include "io/lmptrajec/lmptrajec.h"
#include "io/cpmdinput/cpmdinput.h"

//TODO: check std::ios_base::sync_with_stdio(false)
namespace Vipster{
    const std::map<IOFmt, IO::Plugin const *const> IOPlugins{
            {IOFmt::XYZ, &IO::XYZ},
            {IOFmt::PWI, &IO::PWInput},
            {IOFmt::PWO, &IO::PWOutput},
            {IOFmt::LMP, &IO::LmpInput},
            {IOFmt::DMP, &IO::LmpTrajec},
            {IOFmt::CPI, &IO::CPInput}
    };
    IO::Data readFile(const std::string &fn, IOFmt fmt);
    IO::Data readFile(const std::string &fn, IOFmt fmt, std::string name);
    bool     writeFile(const std::string &fn, IOFmt fmt, const Molecule &m,
                       const BaseParam *p=nullptr,
                       const BaseConfig *c=nullptr,
                       IO::State state={});
}

#endif // IOWRAPPER

