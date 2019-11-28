#ifndef PRESETS_H
#define PRESETS_H

#include "nlohmann/json.hpp"

#include <string>
#include <map>
#include <memory>

namespace Vipster::IO{

class BasePreset
{
public:
    virtual const struct Plugin* getFmt() const = 0;
    virtual std::unique_ptr<BasePreset> copy() const = 0;
    virtual void parseJson(const nlohmann::json&) = 0;
    virtual nlohmann::json toJson() const = 0;
    virtual ~BasePreset() = default;
};

void to_json(nlohmann::json& j, const BasePreset& p);
void from_json(const nlohmann::json& j, BasePreset& p);

constexpr const char* PresetsAbout =
    "IO-presets are used to control HOW the data is "
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

#endif // PRESETS_H
