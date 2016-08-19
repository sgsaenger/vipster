#include "iowrapper.h"

using namespace Vipster;

std::tuple<Molecule,optional<Param>>  Vipster::readFile(std::string fn, IOFmt fmt)
{
    std::ifstream file{fn};
    if(!file){
        throw IOError("Could not open "+fn);
    }
    switch(fmt){
    case IOFmt::xyz:
        return IO::XYZ.parser(fn,file);
        break;
    }
}

void  Vipster::writeFile(const Molecule &m, std::string fn, IOFmt fmt, Param p)
{

}
