#ifndef IOPLUGIN_H
#define IOPLUGIN_H

#include <molecule.h>
#include <fstream>

namespace Vipster{
namespace IO {
    struct BaseParam{
        virtual ~BaseParam() = default;
    };
    struct BaseData{
        Molecule mol{"",0};
        virtual ~BaseData() = default;
    };
}
    struct IOPlugin{
        std::string name;
        std::string extension;
        std::string argument;
        IO::BaseData (*parser)(std::string name, std::ifstream &file);
        bool   (*writer)(const IO::BaseData& d, std::ofstream &file);
    };
    class IOError: public std::runtime_error
    {
        public:
            IOError(std::string reason):std::runtime_error(reason){};
    };
}

#endif // IOPLUGIN_H
