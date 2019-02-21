#ifndef XYZ_CONF_H
#define XYZ_CONF_H

#include "../plugin.h"

namespace Vipster::IO{

struct XYZPreset final: BasePreset{
    enum class Data{None, Charge, Forces};
    enum class Mode{Step, Trajec, Cell};
    Mode filemode;
    Data atomdata;
    XYZPreset(Mode=Mode::Step, Data=Data::None);
    IOFmt getFmt() const override;
    std::unique_ptr<BasePreset> copy() const override;
    void parseJson(const nlohmann::json&) override;
    nlohmann::json toJson() const override;
};

}

#endif // XYZ_CONF_H
