#ifndef IOPLUGIN_H
#define IOPLUGIN_H

#include "molecule.h"
#include <fstream>

namespace Vipster{
namespace IO {
    struct BaseParam{
        std::string name;
        virtual ~BaseParam() = default;
    };
    struct BaseConfig{
        virtual ~BaseConfig() = default;
    };
    struct BaseData{
        Molecule mol{"",0};
        std::unique_ptr<BaseParam> param{};
    };
}
    struct IOPlugin{
        std::string name;
        std::string extension;
        std::string argument;
        IO::BaseData (*parser)(std::string name, std::ifstream &file);
        bool   (*writer)(const Molecule& m, std::ofstream &file,
                         const IO::BaseParam *const p, const IO::BaseConfig *const c);
    };
    class IOError: public std::runtime_error
    {
        public:
            IOError(std::string reason):std::runtime_error(reason){}
    };
}

#endif // IOPLUGIN_H
