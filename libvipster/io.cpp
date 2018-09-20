#include "io.h"
#include <iostream>

using namespace Vipster;

IO::Data Vipster::readFile(const std::string &fn, IOFmt fmt, std::string name)
{
    std::ifstream file{fn};
    try{
        if(!file){
            throw IO::Error("Could not open "+fn);
        }
        if(IOPlugins.find(fmt)==IOPlugins.end()){
            throw IO::Error("Unknown format");
        }
        return IOPlugins.at(fmt)->parser(name, file);
    }
    catch(IO::Error& e){
        std::cout << e.what() << std::endl;
        throw e;
    }
}

IO::Data Vipster::readFile(const std::string &fn, IOFmt fmt)
{
    return readFile(fn, fmt, fn);
}

bool  Vipster::writeFile(const std::string &fn, IOFmt fmt, const Molecule &m,
                         const IO::BaseParam *const p,
                         const IO::BaseConfig *const c,
                         IO::State state)
{
    std::ofstream file{fn};
    if(state.index == -1ul){
        state.index = m.getNstep()-1;
    }
    try{
        if(!file){
            throw IO::Error("Could not open "+fn);
        }
        if(IOPlugins.find(fmt)==IOPlugins.end()){
            throw IO::Error("Unknown format");
        }
        const IO::Plugin * const w = IOPlugins.at(fmt);
        if(w->writer == nullptr){
            throw IO::Error("Read-only format");
        }
        return w->writer(m, file, p, c, state);
    }
    catch(IO::Error &e){
        std::cout << e.what() << std::endl;
        throw e;
    }
}
