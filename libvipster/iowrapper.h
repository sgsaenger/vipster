#ifndef IOWRAPPER
#define IOWRAPPER

#include "definitions.h"
#include "molecule.h"
#include <experimental/optional>
#include <cstdio>
#include "ioplugin.h"
#include "ioplugins/xyz.h"

namespace Vipster{
    std::tuple<Molecule,optional<Param>>  readFile(std::string fn, std::string fmt);
    void        writeFile(const Molecule &m, std::string fn, std::string fmt, Param p);
    const std::map<std::string,IOPlugin> IOPlugins {{IO::XYZ.name,IO::XYZ}};
}

#endif // IOWRAPPER

