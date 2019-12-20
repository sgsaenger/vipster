#ifndef PRESETS_H
#define PRESETS_H

#include "nlohmann/json.hpp"

#include "custommap.h"

#include <string>
#include <variant>

namespace Vipster::IO{

class CustomEnum: public CustomMap<int, std::string>
{
public:
    int value;
    CustomEnum(int value, const std::vector<std::string> &names);
    operator int() const;
    operator const std::string&() const;
    CustomEnum& operator=(int);
    CustomEnum& operator=(const std::string&);
};

using PresetValue = std::variant<bool, uint8_t, CustomEnum>;

class Preset: public CustomMap<std::string, PresetValue>
{
public:
    enum ValIdx { i_bool, i_uint, i_enum };
    const struct Plugin* getFmt() const;
// constructors/destructor
    Preset(const struct Plugin* fmt=nullptr,
               CustomMap<std::string, PresetValue> &&values={});
    Preset(const Preset &) = default;
    Preset(Preset &&) = default;
    Preset& operator=(const Preset &) = default;
    Preset& operator=(Preset &&) = default;
    virtual ~Preset() = default;
private:
    const struct Plugin *fmt;
};

void to_json(nlohmann::json& j, const Preset& p);
void from_json(const nlohmann::json& j, Preset& p);

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
