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
    virtual IOFmt getFmt() const = 0;
    virtual std::unique_ptr<BaseParam> copy() const = 0;
    virtual void parseJson(const nlohmann::json::iterator&) = 0;
    virtual nlohmann::json toJson() const = 0;
    virtual ~BaseParam() = default;
protected:
    BaseParam(std::string);
    BaseParam(const BaseParam&) = default;
    BaseParam(BaseParam&&) = default;
    BaseParam& operator=(const BaseParam&) = default;
    BaseParam& operator=(BaseParam&&) = default;
};

using Parameters = std::map<IOFmt, std::map<std::string, std::unique_ptr<BaseParam>>>;
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

namespace Vipster {

extern IO::Parameters params;

}

#endif // PARAMETERS_H
