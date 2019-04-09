#ifndef GLOBAL_H
#define GLOBAL_H

#include <string>
#include <stdexcept>
#include <array>

namespace Vipster {
    constexpr float pi = 3.14159265358979323846264338327950288f;
    constexpr float rad2deg = 180.f / pi;
    constexpr float deg2rad = pi / 180.f;
    constexpr float bohrrad = 0.52917721092f;
    constexpr float invbohr = 1/bohrrad;

    using ColVec = std::array<uint8_t, 4>;
    using SizeVec = std::array<size_t, 3>;

    class Error:public std::logic_error{
    public:
        Error(std::string reason):std::logic_error(reason){}
    };
}

#endif // GLOBAL_H
