#ifndef CPI_CONF_H
#define CPI_CONF_H

#include "../plugin.h"

namespace Vipster {
namespace IO {

struct CPConfig: BaseConfig{
    enum class Scale{None, Scale, Cartesian};
    bool angstrom;
    Scale scale;
    CPConfig(std::string="", bool=false, Scale=Scale::Scale);
    std::unique_ptr<BaseConfig> copy() override;
    void parseJson(const nlohmann::json&) override;
};

void to_json(nlohmann::json& j,const CPConfig& p);
void from_json(const nlohmann::json& j, CPConfig& p);

}
}

#endif // CPI_CONF_H
