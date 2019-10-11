#include "io.h"
#include <iostream>

using namespace Vipster;

IO::Data Vipster::readFile(const std::string &fn)
{
    // check if file can be read
    std::ifstream file{fn};
    if(!file){
        throw IO::Error("Could not open "+fn);
    }
    // determine filetype
    static const auto IOExt2Fmt = [](){
        std::map<std::string, IOFmt> ext;
        for(const auto& pair: IOPlugins){
            ext[pair.second->extension] = pair.first;
        }
        return ext;
    }();
    std::string name = fn;
    auto pos = fn.find_last_of('.');
    if(pos != fn.npos){
        name = fn.substr(pos+1);
    }else{
        pos = fn.find_last_of("/\\");
        if(pos != fn.npos){
            name = fn.substr(pos+1);
        }
    }
    auto ext = IOExt2Fmt.find(name);
    if(ext == IOExt2Fmt.end()){
        throw IO::Error{"Could not deduce format of file "+fn+
                        "\nPlease specify format explicitely", false};
    }
    // try to parse
    auto tmp = IOPlugins.at(ext->second)->parser(fn, file);
    if(!tmp.mol.getNstep()){
        throw IO::Error("No Molecule could be parsed");
    }
    // return if successful
    return tmp;
}

IO::Data Vipster::readFile(const std::string &fn, IOFmt fmt)
{
    // check if file can be read
    std::ifstream file{fn};
    if(!file){
        throw IO::Error("Could not open "+fn);
    }
    // check filetype
    if(IOPlugins.find(fmt)==IOPlugins.end()){
        throw IO::Error("Unknown format");
    }
    // try to parse
    auto tmp = IOPlugins.at(fmt)->parser(fn, file);
    if(!tmp.mol.getNstep()){
        throw IO::Error("No Molecule could be parsed");
    }
    // return if successful
    return tmp;
}

bool  Vipster::writeFile(const std::string &fn, IOFmt fmt, const Molecule &m,
                         const IO::BaseParam *const p,
                         const IO::BaseConfig *const c,
                         size_t idx)
{
    std::ofstream file{fn};
    if(idx == -1ul){
        idx = m.getNstep()-1;
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
        return w->writer(m, file, p, c, idx);
    }
    catch(IO::Error &e){
        std::cout << e.what() << std::endl;
        throw;
    }
}
