#ifndef IOWRAPPER
#define IOWRAPPER

#include "definitions.h"
#include "molecule.h"
#include <experimental/optional>
#include <fstream>

namespace Vipster{
    using std::experimental::optional;
    typedef std::map<std::string,std::string> Param;
    struct IOPlugin{
        std::string name;
        std::string extension;
        std::string argument;
        std::map<std::string,Param>  *param;
        std::tuple<Molecule,optional<Param>>  (*parser)(std::string fn, std::ifstream &file);
        void        (*writer)(const Molecule &m,std::ofstream &file, Param p);
    };
    class IOError: public std::runtime_error
    {
        public:
            IOError(std::string reason):std::runtime_error(reason){};
    };
}

#include "ioplugins/xyz.h"

namespace Vipster{
    enum class  IOFmt {xyz};

    std::tuple<Molecule,optional<Param>>  readFile(std::string fn, IOFmt fmt);
    void        writeFile(const Molecule &m, std::string fn, IOFmt fmt, Param p);
}

#endif // IOWRAPPER

