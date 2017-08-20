#ifndef IOWRAPPER
#define IOWRAPPER

#include <molecule.h>
#include <ioplugin.h>
#include <map>
#include <ioplugins/xyz.h>
#include <ioplugins/pwinput.h>

namespace Vipster{
    enum class  IOFmt{XYZ,PWI};
    const       std::map<IOFmt, IOPlugin const *const> IOPlugins{
                    {IOFmt::XYZ, &IO::XYZ},
                    {IOFmt::PWI, &IO::PWInput}
                };
    std::shared_ptr<IO::BaseData> readFile(std::string fn, IOFmt fmt);
    std::shared_ptr<IO::BaseData> readFile(std::string fn, IOFmt fmt, std::string name);
    bool        writeFile(std::string fn, IOFmt fmt, const Molecule &m, const IO::BaseParam* const p=nullptr);
}

#endif // IOWRAPPER

