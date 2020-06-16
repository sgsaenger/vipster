#ifndef LMPFF_H
#define LMPFF_H

#include "vipster/molecule.h"
#include "vipster/plugins/lmpinput.h"

struct ForceField{
    std::optional<std::string> pair{};
    std::optional<std::string> bond{};
    std::optional<std::string> angle{};
    std::optional<std::string> dihedral{};
    std::optional<std::string> improper{};
    std::vector<std::string> extra_cmds{};
    std::function<Vipster::Molecule(const Vipster::Step&, const std::string&)> prepareStep{};
    std::function<Vipster::Parameter(const Vipster::Step&)> prepareParameters{};
};

using ForceFields = std::map<std::string, const ForceField*>;
ForceFields defaultForcefields();

#endif // LMPFF_H
