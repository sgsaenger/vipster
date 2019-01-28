#ifndef CONFIGS_H
#define CONFIGS_H

#include "fmt.h"
#include "json.hpp"

#include <string>
#include <map>
#include <memory>

namespace Vipster::IO{

class BaseConfig
{
public:
    std::string name;
    virtual std::unique_ptr<BaseConfig> copy() = 0;
    virtual void parseJson(const nlohmann::json::iterator&) = 0;
    virtual nlohmann::json toJson() = 0;
    virtual ~BaseConfig() = default;
protected:
    BaseConfig(std::string);
    BaseConfig(const BaseConfig&) = default;
    BaseConfig(BaseConfig &&) = default;
    BaseConfig& operator=(const BaseConfig&) = default;
    BaseConfig& operator=(BaseConfig&&) = default;
};

using Configs = std::map<IOFmt, std::map<std::string, std::unique_ptr<BaseConfig>>>;

}

namespace Vipster{
extern IO::Configs configs;
}

#endif // CONFIGS_H
