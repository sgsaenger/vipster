#ifndef IOPLUGIN_H
#define IOPLUGIN_H

#include <molecule.h>
#include <fstream>

#define BUFFLEN 32768

namespace Vipster{
    enum class IOType{ None, Param};
    class IOBase{};
    struct IOPlugin{
        std::string name;
        std::string extension;
        std::string argument;
        std::tuple<Molecule, IOType, IOBase*>  (*parser)(std::string fn, std::ifstream &file);
        void        (*writer)(const Molecule &m,std::ofstream &file, IOBase &p);
    };
    class IOError: public std::runtime_error
    {
        public:
            IOError(std::string reason):std::runtime_error(reason){};
    };
}

#endif // IOPLUGIN_H
