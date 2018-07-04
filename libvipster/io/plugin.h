#ifndef IOPLUGIN_H
#define IOPLUGIN_H

#include "../molecule.h"
#include "../parameters.h"
#include "../configs.h"
#include "fmt.h"

#include <string>
#include <map>
#include <memory>
#include <fstream>

namespace Vipster{

namespace IO {

struct Data{
    Molecule mol{"",0};
    IOFmt fmt;
    std::unique_ptr<BaseParam> param{};
};

struct State{
    size_t index = -1ul;
    AtomFmt atom_fmt = AtomFmt::Crystal;
    CdmFmt cell_fmt = CdmFmt::Bohr;
};

struct Plugin{
    enum Args:uint8_t{None, Param, Config};
    std::string name;
    std::string extension;
    std::string command;
    uint8_t     arguments;
    Data        (*parser)(const std::string& name, std::ifstream &file);
    bool        (*writer)(const Molecule& m, std::ofstream &file,
                          const BaseParam *const p,
                          const BaseConfig *const c,
                          State state);
};

class Error: public std::runtime_error
{
    public:
        Error(const std::string& reason):std::runtime_error(reason){}
};

}
}

#endif // IOPLUGIN_H
