#include "iowrapper.h"
#include <iostream>

using namespace Vipster;

IO::Data Vipster::readFile(std::string fn, IOFmt fmt, std::string name)
{
    std::ifstream file{fn};
    try{
        if(!file){
            throw IO::Error("Could not open "+fn);
        }
        if(IOPlugins.find(fmt)==IOPlugins.end()){
            throw IO::Error("Unknown format");
        }else{
            return IOPlugins.at(fmt)->parser(name, file);
        }
    }
    catch(IO::Error& e){
        std::cout << e.what() << std::endl;
        throw e;
    }
}

IO::Data Vipster::readFile(std::string fn, IOFmt fmt)
{
    return readFile(fn, fmt, fn);
}

bool  Vipster::writeFile(std::string fn, IOFmt fmt, const Molecule &m,
                         const BaseParam *const p,
                         const BaseConfig *const c)
{
    std::ofstream file{fn};
    try{
        if(!file){
            throw IO::Error("Could not open "+fn);
        }
        if(IOPlugins.find(fmt)==IOPlugins.end()){
            throw IO::Error("Unknown format");
        }else{
            const IO::Plugin * const w = IOPlugins.at(fmt);
            if(!w->writer) throw IO::Error("Read-only format");
            return w->writer(m, file, p, c);
        }
    }
    catch(IO::Error &e){
        std::cout << e.what() << std::endl;
        throw e;
    }
}
