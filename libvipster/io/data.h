#ifndef IODATA_H
#define IODATA_H

#include "../molecule.h"
#include "../data.h"
#include "io/parameters.h"

#include <vector>
#include <optional>

namespace Vipster::IO{
    struct Data{
        Molecule mol{"",0};
        std::optional<BaseParam> param{};
        std::vector<std::unique_ptr<const BaseData>> data{};
    };
}

#endif // IODATA_H
