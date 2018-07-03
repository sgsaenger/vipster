#ifndef CPMDINPUT_H
#define CPMDINPUT_H

#include "../ioplugin.h"

namespace Vipster {
namespace IO {

extern const IO::Plugin CPInput;

struct CPParam: BaseParam{
    using Section = std::vector<std::string>;
    Section info;
    Section cpmd;
    Section system;
    Section pimd;
    Section path;
    Section ptddft;
    Section atoms;
    Section dft;
    Section prop;
    Section resp;
    Section linres;
    Section tddft;
    Section hardness;
    Section classic;
    Section vdw;
    Section qmmm;
    CPParam(std::string="", Section={}, Section={}, Section={}, Section={},
            Section={}, Section={}, Section={}, Section={}, Section={}, Section={},
            Section={}, Section={}, Section={}, Section={}, Section={}, Section={});
    std::unique_ptr<BaseParam> copy() override;
};

}
}

#endif // CPMDINPUT_H
