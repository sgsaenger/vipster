#ifndef IOPLUGIN_H
#define IOPLUGIN_H

#include <molecule.h>
#include <fstream>

namespace Vipster{
    enum class IOType{ None, PWParam};
    struct IOBase{};
    struct IOData{
        std::shared_ptr<Molecule> mol;
        IOType data_type;
        std::shared_ptr<IOBase> data;
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
