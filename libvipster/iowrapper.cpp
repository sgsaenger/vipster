#include "iowrapper.h"

using namespace Vipster;

std::tuple<Molecule,optional<Param>>  Vipster::readFile(std::string fn, std::string fmt)
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

void  Vipster::writeFile(const Molecule &m, std::string fn, std::string fmt, Param p)
{

}
