#ifndef IOWRAPPER
#define IOWRAPPER

#include <map>
#include "molecule.h"
#include "ioplugin.h"
#include "ioplugins/xyz.h"
#include "ioplugins/pwinput.h"
#include "ioplugins/pwoutput.h"
#include "ioplugins/lmpinput.h"
#include "ioplugins/lmptrajec.h"
#include "ioplugins/cpmdinput.h"

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

