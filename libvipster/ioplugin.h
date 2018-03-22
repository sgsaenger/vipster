#ifndef IOPLUGIN_H
#define IOPLUGIN_H

#include "molecule.h"
#include <fstream>

namespace Vipster{

enum class IOFmt{XYZ, PWI, PWO, LMP, DMP};

namespace IO {

    struct BaseParam{
        std::string name;
        virtual ~BaseParam() = default;
    };

    struct BaseConfig{
        std::string name;
        virtual ~BaseConfig() = default;
    };

    struct Data{
        Molecule mol{"",0};
        IOFmt fmt;
        std::unique_ptr<BaseParam> param{};
    };
}
    struct IOPlugin{
        enum Args:uint8_t{None, Param, Config};
        std::string name;
        std::string extension;
        std::string command;
        uint8_t     arguments;
        IO::Data    (*parser)(std::string name, std::ifstream &file);
        bool        (*writer)(const Molecule& m, std::ofstream &file,
                              const IO::BaseParam *const p,
                              const IO::BaseConfig *const c);
    };
    class IOError: public std::runtime_error
    {
        public:
            IOError(std::string reason):std::runtime_error(reason){}
    };
}

#endif // IOPLUGIN_H
