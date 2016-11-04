#ifndef PARAM_H
#define PARAM_H

#include <ioplugin.h>

namespace Vipster {
    struct Param: public IOBase
    {
        std::string type;
    };
}

#endif // PARAM_H
