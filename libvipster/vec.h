#ifndef VIPSTER_DEFINITIONS_H
#define VIPSTER_DEFINITIONS_H

#include <array>

namespace Vipster{
    typedef std::array<float,3> Vec;
    bool  operator==(const Vec &v1, const Vec &v2);
    Vec operator+=(Vec &v1, const Vec &v2);
    Vec operator+(Vec v1, const Vec &v2);
    Vec operator-(const Vec &v);
    Vec operator-=(Vec &v1, const Vec &v2);
    Vec operator-(Vec v1, const Vec &v2);
    Vec operator*=(Vec &v, const float &f);
    Vec operator*(Vec v, const float &f);
    Vec operator*(const float &f, Vec v);
    Vec operator/=(Vec &v, const float &f);
    Vec operator/(Vec v, const float &f);
    Vec operator/(const float &f, Vec v);
    float Vec_length(const Vec &v);
    float Vec_dot(const Vec &v1, const Vec &v2);
    Vec Vec_cross(const Vec &v1, const Vec &v2);
}

#endif // VIPSTER_DEFINITIONS_H
