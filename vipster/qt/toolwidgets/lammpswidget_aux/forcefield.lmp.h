#ifndef LMPFF_H
#define LMPFF_H

#include "molecule.h"

struct ForceField{
    std::optional<std::string> required_pair{};
    std::optional<std::string> required_bond{};
    std::optional<std::string> required_angle{};
    std::optional<std::string> required_dihedral{};
    std::optional<std::string> required_improper{};
    std::vector<std::string> extra_cmds{};
    std::function<Vipster::Molecule(const Vipster::Step&, const std::string&)> prepareStep{};
    std::function<void()> prepareParameters{};
};

using ForceFields = std::map<std::string, const ForceField*>;
ForceFields defaultForcefields();

#endif // LMPFF_H
