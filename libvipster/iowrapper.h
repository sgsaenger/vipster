#ifndef IOWRAPPER
#define IOWRAPPER

#include <molecule.h>
#include <ioplugin.h>
//#include <ioplugins/xyz.hpp>
#include <map>

namespace Vipster{
    enum class  IOFmt{XYZ};
    const       std::map<IOFmt,IOPlugin> IOPlugins;
    IOData      readFile(std::string fn, IOFmt fmt);
    bool        writeFile(const IOData &d, std::string fn, IOFmt fmt);
}

#endif // IOWRAPPER

