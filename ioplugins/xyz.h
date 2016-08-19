#ifndef XYZ_H
#define XYZ_H

#include "definitions.h"
#include "iowrapper.h"

namespace Vipster{
    using std::experimental::optional;
    namespace IO{
        std::tuple<Molecule,optional<Param>> XYZ_parser(std::string fn, std::ifstream &file);
        void    XYZ_writer(const Molecule &m, std::ofstream &file,Param p);
        const IOPlugin XYZ{"xyz","xyz","xyz",nullptr,XYZ_parser,XYZ_writer};
    }
}

#endif // XYZ_H
