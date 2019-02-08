#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "fmt.h"
#include "json.hpp"

#include <string>
#include <map>
#include <memory>

namespace Vipster::IO{

class BaseParam
{
public:
    std::string name;
    virtual std::unique_ptr<BaseParam> copy() = 0;
    virtual void parseJson(const nlohmann::json::iterator&) = 0;
    virtual nlohmann::json toJson() = 0;
    virtual ~BaseParam() = default;
protected:
    BaseParam(std::string);
    BaseParam(const BaseParam&) = default;
    BaseParam(BaseParam&&) = default;
    BaseParam& operator=(const BaseParam&) = default;
    BaseParam& operator=(BaseParam&&) = default;
};

using Parameters = std::map<IOFmt, std::map<std::string, std::unique_ptr<BaseParam>>>;

}

namespace Vipster {

extern IO::Parameters params;

}

#endif // PARAMETERS_H
