#ifndef IOPLUGIN_H
#define IOPLUGIN_H

#include "../molecule.h"
#include "../settings.h"
#include "parameters.h"
#include "presets.h"
#include "data.h"

#include <string>
#include <map>
#include <memory>
#include <fstream>
#include <iomanip>

namespace Vipster::IO{

struct Plugin{
    // TODO: remove Args?
    enum Args:uint8_t{None=0x0, Read=0x1, Write=0x2, Param=0x4, Preset=0x8};
    std::string name;
    std::string extension;
    std::string command;
    uint8_t     arguments;
    std::function<IO::Data(const std::string& name, std::istream &file)> parser{};
    std::function<bool(const Molecule& m, std::ostream &file,
                       const BaseParam *const p, const BasePreset *const c,
                       size_t idx)> writer{};
    std::function<std::unique_ptr<BaseParam>()> makeParam{};
    std::function<std::unique_ptr<BasePreset>()> makePreset{};
};

class Error: public std::runtime_error
{
    public:
        Error(const std::string& reason, bool fatal=true)
            : std::runtime_error{reason}, fatal{fatal}
        {}
        bool fatal;
};

}

#endif // IOPLUGIN_H
