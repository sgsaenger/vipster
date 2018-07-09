#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "fmt.h"
#include "json.hpp"

#include <string>
#include <map>
#include <memory>

namespace Vipster {
namespace IO {

class BaseParam
{
public:
    std::string name;
    virtual std::unique_ptr<BaseParam> copy() = 0;
    virtual void parseJson(const nlohmann::json&) = 0;
    virtual nlohmann::json toJson() = 0;
    virtual ~BaseParam() = default;
protected:
    BaseParam(std::string);
};

using Parameters = std::multimap<IOFmt, std::unique_ptr<BaseParam>>;

}

extern IO::Parameters params;

}

#endif // PARAMETERS_H
