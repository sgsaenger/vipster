#include "iowrapper.h"
#include <iostream>

using namespace Vipster;

IO::BaseData Vipster::readFile(std::string fn, IOFmt fmt, std::string name)
{
    std::ifstream file{fn};
    try{
        if(!file){
            throw IOError("Could not open "+fn);
        }
        if(IOPlugins.find(fmt)==IOPlugins.end()){
            throw IOError("Unknown format");
        }else{
            return IOPlugins.at(fmt)->parser(name, file);
        }
    }
    catch(IOError& e){
        std::cout << e.what() << std::endl;
        throw e;
    }
}

IO::BaseData Vipster::readFile(std::string fn, IOFmt fmt)
{
    return readFile(fn, fmt, fn);
}

bool  Vipster::writeFile(std::string fn, IOFmt fmt, const Molecule &m,
                         const IO::BaseParam *const p,
                         const IO::BaseConfig *const c)
{
    std::ofstream file{fn};
    try{
        if(!file){
            throw IOError("Could not open "+fn);
        }
        if(IOPlugins.find(fmt)==IOPlugins.end()){
            throw IOError("Unknown format");
        }else{
            const IOPlugin * const w = IOPlugins.at(fmt);
            if(!w->writer) throw IOError("Read-only format");
            return w->writer(m, file, p, c);
        }
    }
    catch(IOError &e){
        std::cout << e.what() << std::endl;
        throw e;
    }
}
