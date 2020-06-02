#ifndef PLUG_PY_H
#define PLUG_PY_H

#include <filesystem>
#include "../global.py.h"
#include "plugin.h"

namespace Vipster::Py{
void Plugins(py::module& io);

class Plugin : public IO::Plugin
{
public:
    static Plugin* create(std::filesystem::path file);
private:
    Plugin()=default;
    IO::Data parser_impl(const std::string& n, std::istream &file);
    bool writer_impl(const Molecule &m, std::ostream &file,
                     const std::optional<IO::Parameter>& p,
                     const std::optional<IO::Preset>& c,
                     size_t idx);
    IO::Parameter makeParam_impl();
    IO::Preset makePreset_impl();
    py::object pyReader;
    py::object pyWriter;
    IO::Parameter::BaseMap param;
    IO::Preset::BaseMap preset;
};
}

#endif // PLUG_PY_H
