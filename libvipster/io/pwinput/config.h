#ifndef PWI_CONF_H
#define PWI_CONF_H

#include "../plugin.h"

namespace Vipster{
namespace IO{

struct PWConfig: BaseConfig{
    enum class AtomFmt {Bohr, Angstrom, Crystal, Alat, Current};
    enum class CellFmt {Bohr, Angstrom, Current};
    AtomFmt atoms;
    CellFmt cell;
    PWConfig(std::string="", AtomFmt=AtomFmt::Current, CellFmt=CellFmt::Current);
    std::unique_ptr<BaseConfig> copy() override;
    void parseJson(const nlohmann::json&) override;
};

void to_json(nlohmann::json& j,const PWConfig& p);
void from_json(const nlohmann::json& j, PWConfig& p);

}
}

#endif // PWI_CONF_H
