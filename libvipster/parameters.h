#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "io/fmt.h"

#include <string>
#include <map>
#include <memory>

namespace Vipster {

class BaseParam
{
public:
    std::string name;
    virtual std::unique_ptr<BaseParam> copy() = 0;
    virtual ~BaseParam() = default;
protected:
    BaseParam(std::string);
};

using Parameters = std::multimap<IOFmt, std::unique_ptr<BaseParam>>;

extern Parameters params;

}

#endif // PARAMETERS_H
