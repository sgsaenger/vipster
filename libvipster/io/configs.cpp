#include "configs.h"
#include "../io.h"

using namespace Vipster;

IO::BaseConfig::BaseConfig(std::string name)
    :name{name}
{}

IO::Configs Vipster::configs = [](){
    IO::Configs tmp;
    for(const auto& p: IOPlugins){
        if(p.second->arguments & IO::Plugin::Args::Config){
            tmp[p.first]["default"] = p.second->makeConfig("default");
        }
    }
    return tmp;
}();
