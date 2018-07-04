#ifndef CONFIGS_H
#define CONFIGS_H

#include "io/fmt.h"

#include <string>
#include <map>
#include <memory>

namespace Vipster {

class BaseConfig
{
public:
    std::string name;
    virtual std::unique_ptr<BaseConfig> copy() = 0;
    virtual ~BaseConfig() = default;
protected:
    BaseConfig(std::string);
};

using Configs = std::multimap<IOFmt, std::unique_ptr<BaseConfig>>;

extern Configs configs;

}

#endif // CONFIGS_H
