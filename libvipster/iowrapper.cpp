#include "iowrapper.h"
#include <boost/spirit/include/qi.hpp>

using namespace Vipster;

IOData Vipster::readFile(std::string fn, IOFmt fmt)
{
    std::ifstream file{fn};
    if(!file){
        throw IOError("Could not open "+fn);
    }
    if(Vipster::IOPlugins.find(fmt)==Vipster::IOPlugins.end()){
        throw IOError("Unknown format");
    }else{
        return Vipster::IOPlugins.at(fmt)->parser(fn,file);
    }
}

//bool  Vipster::writeFile(const IOData &d, std::string fn, IOFmt fmt)
//{
//    return false;
//}
