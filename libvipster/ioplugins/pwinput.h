#ifndef PWINPUT_H
#define PWINPUT_H

#include "ioplugin.h"

namespace Vipster{
namespace IO{
extern const IOPlugin PWInput;
struct PWParam: BaseParam{
    std::map<std::string, std::string> control;
    std::map<std::string, std::string> system;
    std::map<std::string, std::string> electrons;
    std::map<std::string, std::string> ions;
    std::map<std::string, std::string> cell;
};
struct PWData: Vipster::IO::BaseData{
    PWParam data;
};
}
}


#endif // PWINPUT_H
