#ifndef IOPLUGIN_H
#define IOPLUGIN_H

#include <molecule.h>
#include <fstream>

namespace Vipster{
namespace IO {
    struct BaseParam{
//        virtual Vipster::IOType getIOType()=0;
        virtual ~BaseParam(){};
    };
}
    struct IOData{
        Molecule mol;
        std::shared_ptr<IO::BaseParam> data;
    };
    struct IOPlugin{
        std::string name;
        std::string extension;
        std::string argument;
        IOData (*parser)(std::string fn, std::ifstream &file);
        bool   (*writer)(const IOData& d, std::ofstream &file);
    };
    class IOError: public std::runtime_error
    {
        public:
            IOError(std::string reason):std::runtime_error(reason){};
    };
}

#endif // IOPLUGIN_H
