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

//TODO: check std::ios_base::sync_with_stdio(false)
namespace Vipster{
    enum class  IOFmt{XYZ, PWI, PWO, LMP, DMP};
    const       std::map<IOFmt, IOPlugin const *const> IOPlugins{
                    {IOFmt::XYZ, &IO::XYZ},
                    {IOFmt::PWI, &IO::PWInput},
                    {IOFmt::PWO, &IO::PWOutput},
                    {IOFmt::LMP, &IO::LmpInput},
                    {IOFmt::DMP, &IO::LmpTrajec}
                };
    IO::BaseData readFile(std::string fn, IOFmt fmt);
    IO::BaseData readFile(std::string fn, IOFmt fmt, std::string name);
    bool        writeFile(std::string fn, IOFmt fmt, const Molecule &m,
                          const IO::BaseParam *const p=nullptr,
                          const IO::BaseConfig *const c=nullptr);
}

#endif // IOWRAPPER

