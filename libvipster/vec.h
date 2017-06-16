#ifndef VIPSTER_DEFINITIONS_H
#define VIPSTER_DEFINITIONS_H

#include <array>
#include <cmath>
#include <limits>
#include <ostream>

namespace Vipster{
typedef std::array<float,3> Vec;
typedef std::array<Vec,3> Mat;

inline std::ostream& operator<<(std::ostream &s, const Vec &v)
{
    s << "Vec[" << v[0] << ", " << v[1] << ", " << v[2] << "]";
    return s;
}

inline bool  operator==(const Vec &v1, const Vec &v2)
{
    return std::abs(v1[0]-v2[0]) < std::numeric_limits<float>::epsilon()&&
           std::abs(v1[1]-v2[1]) < std::numeric_limits<float>::epsilon()&&
           std::abs(v1[2]-v2[2]) < std::numeric_limits<float>::epsilon();
}

inline bool  operator!=(const Vec &v1, const Vec &v2)
{
    return !(v1==v2);
}

inline Vec operator +=(Vec &v1, const Vec &v2)
{
    v1[0]+=v2[0];
    v1[1]+=v2[1];
    v1[2]+=v2[2];
    return v1;
}

inline Vec operator +(Vec v1, const Vec &v2)
{
    return v1+=v2;
}

inline Vec operator +=(Vec &v, const float &f)
{
    v[0]+=f;
    v[1]+=f;
    v[2]+=f;
    return v;
}

inline Vec operator +(Vec v, const float &f)
{
    return v+=f;
}

inline Vec operator +(const float &f, Vec v)
{
    return v+=f;
}

inline Vec operator -(const Vec &v)
{
    return {{-v[0],-v[1],-v[2]}};
}

inline Vec operator -=(Vec &v1, const Vec &v2)
{
    v1[0]-=v2[0];
    v1[1]-=v2[1];
    v1[2]-=v2[2];
    return v1;
}

inline Vec operator -(Vec v1, const Vec &v2)
{
    return v1-=v2;
}

inline Vec operator -=(Vec &v, const float &f)
{
    v[0]-=f;
    v[1]-=f;
    v[2]-=f;
    return v;
}

inline Vec operator -(Vec v, const float &f)
{
    return v-=f;
}

inline Vec operator *=(Vec &v, const float &f)
{
    v[0]*=f;
    v[1]*=f;
    v[2]*=f;
    return v;
}

inline Vec operator *(Vec v, const float &f)
{
    return v*=f;
}

inline Vec operator *(const float &f, Vec v)
{
    return v*=f;
}

inline Vec operator /=(Vec &v, const float &f)
{
    v[0]/=f;
    v[1]/=f;
    v[2]/=f;
    return v;
}

inline Vec operator /(Vec v, const float &f)
{
    return v/=f;
}

inline float Vec_dot(const Vec &v1, const Vec &v2)
{
    return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
}

inline float Vec_length(const Vec &v)
{
    return std::sqrt(Vec_dot(v,v));
}

inline Vec Vec_cross(const Vec &v1, const Vec &v2)
{
    return {{v1[1]*v2[2]-v1[2]*v2[1],
             v1[2]*v2[0]-v1[0]*v2[2],
             v1[0]*v2[1]-v1[1]*v2[0]}};
}

inline Vec operator*(const Mat &m, const Vec &v)
{
    Vec t;
    t[0] = Vec_dot(m[0],v);
    t[1] = Vec_dot(m[1],v);
    t[2] = Vec_dot(m[2],v);
    return t;
}

inline Mat Mat_trans(const Mat &m)
{
    return {{Vec{{m[0][0],m[1][0],m[2][0]}},
             Vec{{m[0][1],m[1][1],m[2][1]}},
             Vec{{m[0][2],m[1][2],m[2][2]}}}};
}

inline Vec operator*(const Vec &v, const Mat&m)
{
    Vec t;
    Mat mt = Mat_trans(m);
    t[0] = Vec_dot(mt[0],v);
    t[1] = Vec_dot(mt[1],v);
    t[2] = Vec_dot(mt[2],v);
    return t;
}

inline float Mat_det(const Mat &m)
{
    return  m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])
           +m[0][1]*(m[1][2]*m[2][0]-m[1][0]*m[2][2])
           +m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]);
}

inline Mat Mat_inv(const Mat &m)
{
    float d = Mat_det(m);
    if(std::abs(d) < std::numeric_limits<float>::epsilon())
    {
        throw std::invalid_argument("Mat_inv: singular matrix has no inverse!");
    }
    d = 1/d;
    Mat inv;
    inv[0][0] = (m[1][1]*m[2][2]-m[2][1]*m[1][2])*d;
    inv[1][0] =-(m[1][0]*m[2][2]-m[2][0]*m[1][2])*d;
    inv[2][0] = (m[1][0]*m[2][1]-m[2][0]*m[1][1])*d;
    inv[0][1] =-(m[0][1]*m[2][2]-m[2][1]*m[0][2])*d;
    inv[1][1] = (m[0][0]*m[2][2]-m[2][0]*m[0][2])*d;
    inv[2][1] =-(m[0][0]*m[2][1]-m[2][0]*m[0][1])*d;
    inv[0][2] = (m[0][1]*m[1][2]-m[1][1]*m[0][2])*d;
    inv[1][2] =-(m[0][0]*m[1][2]-m[1][0]*m[0][2])*d;
    inv[2][2] = (m[0][0]*m[1][1]-m[1][0]*m[0][1])*d;
    return inv;
}

inline std::ostream& operator<<(std::ostream &s, const Mat &v)
{
    s << "Mat[[" << v[0][0] << ", " << v[0][1] << ", " << v[0][2] << "]\n"
      << "    [" << v[1][0] << ", " << v[1][1] << ", " << v[1][2] << "]\n"
      << "    [" << v[2][0] << ", " << v[2][1] << ", " << v[2][2] << "]]";
    return s;
}

}

#endif // VIPSTER_DEFINITIONS_H