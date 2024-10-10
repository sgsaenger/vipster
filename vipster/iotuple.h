#ifndef IODATA_H
#define IODATA_H

#include "molecule.h"
#include "data.h"
#include "parameters.h"

#include <vector>
#include <optional>

namespace Vipster{
    // Simple tuple of Types produced from parsing files
    // Implemented as a class because MSVC is buggy...
    // for some reason, it insists on using the deleted copy constructor in vector.push_back
    class IOTuple: public std::tuple<Molecule,
                               std::optional<Parameter>,
                               DataList>
    {
    public:
        using std::tuple<Molecule, std::optional<Parameter>, DataList>::tuple;
        IOTuple(const IOTuple&) = delete;
        IOTuple& operator=(const IOTuple&) = delete;
        IOTuple(IOTuple&&) = default;
        IOTuple& operator=(IOTuple&&) = default;
    };
}

namespace std
{
    template<>
    struct tuple_size<Vipster::IOTuple>: std::integral_constant<size_t, 3>{};
    template<>
    struct tuple_element<0, Vipster::IOTuple> { using type = Vipster::Molecule; };
    template<>
    struct tuple_element<1, Vipster::IOTuple> { using type = std::optional<Vipster::Parameter>; };
    template<>
    struct tuple_element<2, Vipster::IOTuple> { using type = Vipster::DataList; };
}

#endif // IODATA_H
