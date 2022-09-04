#ifndef IOPLUGIN_H
#define IOPLUGIN_H

#include "molecule.h"
#include "settings.h"
#include "parameters.h"
#include "presets.h"
#include "iotuple.h"

#include <string>
#include <map>
#include <memory>
#include <fstream>
#include <iomanip>

namespace Vipster{

using IOTuple = std::tuple<Molecule,
                           std::optional<Parameter>,
                           DataList>;

struct Plugin{
    std::string name;
    std::string extension;
    std::string command;
    std::function<IOTuple(const std::string& name, std::istream &file)> parser{};
    std::function<bool(const Molecule& m, std::ostream &file,
                       const std::optional<Parameter>& p,
                       const std::optional<Preset>& c,
                       size_t idx)> writer{};
    std::function<Parameter()> makeParam{};
    std::function<Preset()> makePreset{};
};

using PluginList = std::vector<const Plugin*>;
PluginList defaultPlugins();

class IOError: public std::runtime_error
{
    public:
        IOError(const std::string& reason)
            : std::runtime_error{reason}
        {}
};

}

#endif // IOPLUGIN_H
