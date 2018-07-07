#ifndef CONFIGS_H
#define CONFIGS_H

#include "io/fmt.h"
#include "json.hpp"

#include <string>
#include <map>
#include <memory>

namespace Vipster {

class BaseConfig
{
public:
    std::string name;
    virtual std::unique_ptr<BaseConfig> copy() = 0;
    virtual void parseJson(const nlohmann::json&) = 0;
    virtual nlohmann::json toJson() = 0;
    virtual ~BaseConfig() = default;
protected:
    BaseConfig(std::string);
};

using Configs = std::multimap<IOFmt, std::unique_ptr<BaseConfig>>;

extern Configs configs;

}

#endif // CONFIGS_H
