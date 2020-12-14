#ifndef PLUG_PY_H
#define PLUG_PY_H

#include <filesystem>
#include "global.py.h"
#include "plugin.h"

namespace Vipster{
std::ostream& operator<<(std::ostream &os, const Plugin *p);
}

namespace Vipster::Py{
void Plugins(py::module& m);

class Plugin : public Vipster::Plugin
{
public:
    static Plugin* create(std::filesystem::path file);
private:
    Plugin()=default;
    IOTuple parser_impl(const std::string& n, std::istream &file);
    bool writer_impl(const Molecule &m, std::ostream &file,
                     const std::optional<Parameter>& p,
                     const std::optional<Preset>& c,
                     size_t idx);
    Parameter makeParam_impl();
    Preset makePreset_impl();
    py::object pyReader;
    py::object pyWriter;
    Parameter::BaseMap param;
    Preset::BaseMap preset;
};
}

#endif // PLUG_PY_H
