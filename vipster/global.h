#ifndef GLOBAL_H
#define GLOBAL_H

#include <string>
#include <stdexcept>
#include <array>

namespace Vipster {
    constexpr double pi = 3.14159265358979323846264338327950288;
    constexpr double rad2deg = 180. / pi;
    constexpr double deg2rad = pi / 180.;
    constexpr double bohrrad = 0.52917721092;
    constexpr double invbohr = 1./bohrrad;

    using ColVec = std::array<uint8_t, 4>;
    using SizeVec = std::array<size_t, 3>;
    using DiffVec = std::array<int, 3>;

    constexpr static std::array<Vipster::ColVec, 5> defaultColors{
        Vipster::ColVec{80, 0, 0, 200},
        Vipster::ColVec{0, 80, 0, 200},
        Vipster::ColVec{80, 80, 0, 200},
        Vipster::ColVec{80, 0, 80, 200},
        Vipster::ColVec{0, 80, 80, 200}
    };

    class Error:public std::logic_error{
    public:
        Error(const std::string &reason):std::logic_error(reason){}
    };
}

#endif // GLOBAL_H
