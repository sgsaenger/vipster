#ifndef IOWRAPPER
#define IOWRAPPER

#include <molecule.h>
#include <ioplugin.h>
#include <ioplugins/xyz.hpp>

namespace Vipster{
    enum class IOFmt{XYZ};
    const std::map<IOFmt,IOPlugin> IOPlugins {{IOFmt::XYZ, IO::XYZ}};
    std::tuple<Molecule, IOType, IOBase*>  readFile(std::string fn, IOFmt fmt);
    void        writeFile(const Molecule &m, std::string fn, IOFmt fmt, IOBase p);
//    const std::map<std::string,IOPlugin> IOPlugins {{IO::XYZ.name,IO::XYZ}};
}

#endif // IOWRAPPER

