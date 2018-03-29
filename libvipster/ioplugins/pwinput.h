#ifndef PWINPUT_H
#define PWINPUT_H

#include "../ioplugin.h"

namespace Vipster{
namespace IO{

extern const IO::Plugin PWInput;

using PWNamelist = std::map<std::string, std::string>;

struct PWParam: BaseParam{
    PWNamelist control;
    PWNamelist system;
    PWNamelist electrons;
    PWNamelist ions;
    PWNamelist cell;
    std::unique_ptr<BaseParam> copy();
};

void to_json(nlohmann::json& j,const PWParam& p);
void from_json(const nlohmann::json& j, PWParam& p);

//TODO: PWConfig? set format-preferences?

}
}


#endif // PWINPUT_H
