#include "iowrapper.h"

using namespace Vipster;

std::tuple<Molecule, IOType, IOBase*>  Vipster::readFile(std::string fn, IOFmt fmt)
{
    std::ifstream file{fn};
    if(!file){
        throw IOError("Could not open "+fn);
    }
    if(Vipster::IOPlugins.find(fmt)==Vipster::IOPlugins.end()){
        throw IOError("Unknown format");
    }else{
        return Vipster::IOPlugins.at(fmt).parser(fn,file);
    }
}

void  Vipster::writeFile(const Molecule &m, std::string fn, IOFmt fmt, IOBase p)
{

}
