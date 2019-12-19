#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "nlohmann/json.hpp"

#include "custommap.h"

#include <string>
#include <variant>

namespace Vipster::IO{

using ParamValue = std::variant<std::string,
                                std::vector<std::string>,
                                std::map<std::string, std::string>>;
class BaseParam : public CustomMap<std::string, ParamValue>
{
public:
    enum ValIdx { i_str, i_strvec, i_strmap };
    const struct Plugin* getFmt() const;
// constructors/destructor
    BaseParam(const struct Plugin* fmt=nullptr,
              CustomMap<std::string, ParamValue> &&values={});
    BaseParam(const BaseParam &) = default;
    BaseParam(BaseParam &&) = default;
    BaseParam& operator=(const BaseParam &) = default;
    BaseParam& operator=(BaseParam &&) = default;
    virtual ~BaseParam() = default;
private:
    const struct Plugin *fmt;
};

void to_json(nlohmann::json& j, const BaseParam& p);
void from_json(const nlohmann::json& j, BaseParam& p);

constexpr const char* ParametersAbout =
    "Parameter presets are used to control additional data that is written "
    "to the formatted target file.\n\n"
    "E.g. PWScf input files contain some Fortran-style namelists that "
    "control the calculation details.\n"
    "While information that is tied to the structure (e.g. number of atoms or "
    "cell parameters) will be extracted upon parsing and reinserted in the correct "
    "format upon writing, decoupled settings (e.g. the type of calculation or "
    "convergence criteria) are saved in a parameter preset to be reused at will.\n"
    "The main use is in automatic conversion on the command-line:\n\n"
    "$vipster convert xyz input.xyz pwi output.pwi -p default"
    ;

}

#endif // PARAMETERS_H
