#ifndef IOPLUGIN_H
#define IOPLUGIN_H

#include "definitions.h"
#include "molecule.h"
#include <experimental/optional>
#include <fstream>
#include <cstring>

#define BUFFLEN 32768

namespace Vipster{
    using std::experimental::optional;
    struct Param{
        std::string type;
    };
    struct IOPlugin{
        std::string name;
        std::string extension;
        std::string argument;
        std::map<std::string,Param>  *param;
        std::tuple<Molecule,optional<Param>>  (*parser)(std::string fn, std::ifstream &file);
        void        (*writer)(const Molecule &m,std::ofstream &file, Param &p);
    };
    class IOError: public std::runtime_error
    {
        public:
            IOError(std::string reason):std::runtime_error(reason){};
    };
}

#endif // IOPLUGIN_H
