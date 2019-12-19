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
    std::string name;
    std::string extension;
    std::string command;
    std::function<IO::Data(const std::string& name, std::istream &file)> parser{};
    std::function<bool(const Molecule& m, std::ostream &file,
                       const std::optional<Parameter>& p,
                       const std::optional<Preset>& c,
                       size_t idx)> writer{};
    std::function<Parameter()> makeParam{};
    std::function<Preset()> makePreset{};
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
