#include <vec.h>
#include <cmath>
#include <limits>

using namespace Vipster;

Vec Vipster::operator +=(Vec &v1, const Vec &v2)
{
    v1[0]+=v2[0];
    v1[1]+=v2[1];
    v1[2]+=v2[2];
    return v1;
}

Vec Vipster::operator +(Vec v1, const Vec &v2)
{
    return v1+=v2;
}

Vec Vipster::operator +=(Vec &v, const float &f)
{
    v[0]+=f;
    v[1]+=f;
    v[2]+=f;
    return v;
}

Vec Vipster::operator +(Vec v, const float &f)
{
    return v+=f;
}

Vec Vipster::operator +(const float &f, Vec v)
{
    return v+=f;
}

Vec Vipster::operator-=(Vec &v1, const Vec &v2)
{
    v1[0]-=v2[0];
    v1[1]-=v2[1];
    v1[2]-=v2[2];
    return v1;
}

Vec Vipster::operator-(Vec v1, const Vec &v2)
{
    return v1-=v2;
}

Vec Vipster::operator-(const Vec &v)
{
    return Vec{{-v[0],-v[1],-v[2]}};
}

Vec Vipster::operator -=(Vec &v, const float &f)
{
    v[0]-=f;
    v[1]-=f;
    v[2]-=f;
    return v;
}

Vec Vipster::operator -(Vec v, const float &f)
{
    return v-=f;
}

Vec Vipster::operator -(const float &f, Vec v)
{
    return v-=f;
}

Vec Vipster::operator*=(Vec &v, const float &f)
{
    v[0]*=f;
    v[1]*=f;
    v[2]*=f;
    return v;
}

Vec Vipster::operator*(Vec v, const float &f)
{
    return v*=f;
}

Vec Vipster::operator*(const float &f, Vec v)
{
    return v*=f;
}

Vec Vipster::operator/=(Vec &v, const float &f)
{
    v[0]/=f;
    v[1]/=f;
    v[2]/=f;
    return v;
}

Vec Vipster::operator/(Vec v, const float &f)
{
    return v/=f;
}

Vec Vipster::operator/(const float &f, Vec v)
{
    return v/=f;
}

bool Vipster::operator ==(const Vec &v1, const Vec &v2)
{
    return (v1[0]-v2[0]) < std::numeric_limits<float>::epsilon()&&
           (v1[1]-v2[1]) < std::numeric_limits<float>::epsilon()&&
           (v1[2]-v2[2]) < std::numeric_limits<float>::epsilon();
}

float Vipster::Vec_length(const Vec &v)
{
    return std::sqrt(Vec_dot(v,v));
}

float Vipster::Vec_dot(const Vec &v1, const Vec &v2)
{
    return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
}

Vec Vipster::Vec_cross(const Vec &v1, const Vec &v2)
{
    return Vec{{v1[1]*v2[2] - v1[2]*v2[1],
                v1[2]*v2[0] - v1[0]*v2[2],
                v1[0]*v2[1] - v1[1]*v2[0]}};
}
