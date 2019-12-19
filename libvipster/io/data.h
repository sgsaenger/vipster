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
        std::optional<Parameter> param{};
        std::vector<std::unique_ptr<const BaseData>> data{};
    // prohibit copy, still needed for pybind11
        Data(const Data&)=delete;
        Data(Data&&)=default;
        Data& operator=(const Data&)=delete;
        Data& operator=(Data&&)=default;
    };
}

#endif // IODATA_H
