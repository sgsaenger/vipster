#include "iowrapper.h"
#include <iostream>
//#include <boost/spirit/include/qi.hpp>

using namespace Vipster;

IO::BaseData Vipster::readFile(std::string fn, IOFmt fmt, std::string name)
{
    std::ifstream file{fn};
    try{
        if(!file){
            throw IOError("Could not open "+fn);
        }
        if(Vipster::IOPlugins.find(fmt)==Vipster::IOPlugins.end()){
            throw IOError("Unknown format");
        }else{
            return Vipster::IOPlugins.at(fmt)->parser(name, file);
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

//bool  Vipster::writeFile(const IO::BaseData &d, std::string fn, IOFmt fmt)
//{
//    return false;
//}
