#ifndef GLOBAL_H
#define GLOBAL_H

#include <stdexcept>

namespace Vipster {
    constexpr float bohrrad = 0.52917721092;
    constexpr float invbohr = 1/bohrrad;

    class Error:public std::logic_error{
    public:
        Error(std::string reason):std::logic_error(reason){}
    };
}

#endif // GLOBAL_H
