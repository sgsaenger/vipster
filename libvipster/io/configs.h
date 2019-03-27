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
    virtual IOFmt getFmt() const = 0;
    virtual std::unique_ptr<BaseConfig> copy() const = 0;
    virtual void parseJson(const nlohmann::json::iterator&) = 0;
    virtual nlohmann::json toJson() const = 0;
    virtual ~BaseConfig() = default;
protected:
    BaseConfig(std::string);
    BaseConfig(const BaseConfig&) = default;
    BaseConfig(BaseConfig &&) = default;
    BaseConfig& operator=(const BaseConfig&) = default;
    BaseConfig& operator=(BaseConfig&&) = default;
};

using Configs = std::map<IOFmt, std::map<std::string, std::unique_ptr<BaseConfig>>>;
constexpr const char* ConfigsAbout =
    "IO-Config presets are used to control HOW the data is "
    "written to the formatted target file.\n\n"
    "E.g. XYZ canonically contains one or more steps of a trajectory, "
    "containing solely the atom types and coordinates.\n"
    "With a certain preset, you can choose whether you want to write "
    "only the step you're working on, or the complete trajectory.\n"
    "Also, you can enable non-standard additional data like atom-charges or force-vectors.\n"
    "The main use is in automatic conversion on the command-line:\n\n"
    "$vipster convert xyz input.xyz xyz output.xyz -c default"
    ;

}

namespace Vipster{
extern IO::Configs configs;
}

#endif // CONFIGS_H
