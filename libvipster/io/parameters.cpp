#include "parameters.h"
#include "../io.h"

using namespace Vipster;

IO::BaseParam::BaseParam(std::string name)
    : name{name}
{}

IO::Parameters Vipster::params = [](){
    IO::Parameters tmp;
    for(const auto& p: IOPlugins){
        if(p.second->arguments & IO::Plugin::Args::Param){
            tmp[p.first]["default"] = p.second->makeParam("default");
        }
    }
    return tmp;
}();
