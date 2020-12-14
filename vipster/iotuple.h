#ifndef IODATA_H
#define IODATA_H

#include "molecule.h"
#include "data.h"
#include "parameters.h"

#include <vector>
#include <optional>

namespace Vipster{
    using IOTuple = std::tuple<Molecule,
                               std::optional<Parameter>,
                               DataList>;
}

#endif // IODATA_H
