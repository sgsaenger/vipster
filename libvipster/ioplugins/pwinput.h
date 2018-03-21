#ifndef PWINPUT_H
#define PWINPUT_H

#include "../ioplugin.h"

namespace Vipster{
namespace IO{

extern const IOPlugin PWInput;

using PWNamelist = std::map<std::string, std::string>;

struct PWParam: BaseParam{
    PWNamelist control;
    PWNamelist system;
    PWNamelist electrons;
    PWNamelist ions;
    PWNamelist cell;
};

}
}


#endif // PWINPUT_H
