#ifndef IOPLUGIN_H
#define IOPLUGIN_H

#include "../molecule.h"
#include "../data.h"
#include "parameters.h"
#include "configs.h"
#include "settings.h"
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

struct Plugin{
    enum Args:uint8_t{None, Param, Config};
    std::string name;
    std::string extension;
    std::string command;
    uint8_t     arguments;
    Data        (*parser)(const std::string& name, std::istream &file);
    bool        (*writer)(const Molecule& m, std::ostream &file,
                          const BaseParam *const p,
                          const BaseConfig *const c,
                          size_t idx) = nullptr;
    std::unique_ptr<BaseParam> (*makeParam)(const std::string& name) = nullptr;
    std::unique_ptr<BaseConfig> (*makeConfig)(const std::string& name) = nullptr;
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
