#ifndef IOPLUGIN_H
#define IOPLUGIN_H

#include "molecule.h"
#include <fstream>

namespace Vipster{
namespace IO {
    constexpr int linelen = 256;
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
        std::shared_ptr<IO::BaseData> (*parser)(std::string name, std::ifstream &file);
        bool   (*writer)(const Molecule& m, std::ofstream &file, const IO::BaseParam* p);
    };
    class IOError: public std::runtime_error
    {
        public:
            IOError(std::string reason):std::runtime_error(reason){}
    };
}

#endif // IOPLUGIN_H
