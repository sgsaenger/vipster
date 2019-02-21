#ifndef IOPLUGIN_H
#define IOPLUGIN_H

#include "../molecule.h"
#include "../data.h"
#include "parameters.h"
#include "presets.h"
#include "fmt.h"

#include <string>
#include <map>
#include <memory>
#include <fstream>
#include <iomanip>

namespace Vipster::IO{

struct Data{
    Molecule mol{"",0};
    std::unique_ptr<BaseParam> param{};
    std::vector<std::unique_ptr<const BaseData>> data{};
};

struct State{
    size_t index = -1ul;
    AtomFmt atom_fmt = AtomFmt::Crystal;
    CdmFmt cell_fmt = CdmFmt::Bohr;
};

struct Plugin{
    enum Args:uint8_t{None=0x0, Param=0x1, Preset=0x2};
    std::string name;
    std::string extension;
    std::string command;
    uint8_t     arguments;
    Data        (*parser)(const std::string& name, std::ifstream &file);
    bool        (*writer)(const Molecule& m, std::ofstream &file,
                          const BaseParam *const p,
                          const BasePreset *const c,
                          State state) = nullptr;
    std::unique_ptr<BaseParam> (*makeParam)() = nullptr;
    std::unique_ptr<BasePreset> (*makePreset)() = nullptr;
};

class Error: public std::runtime_error
{
    public:
        Error(const std::string& reason):std::runtime_error(reason){}
};

}

#endif // IOPLUGIN_H
