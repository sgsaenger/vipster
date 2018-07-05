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
};

const CPConfig CPConfigDefault{
    "default",
    false,
    CPConfig::Scale::None
};

}
}

#endif // CPI_CONF_H
