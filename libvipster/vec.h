#ifndef VIPSTER_DEFINITIONS_H
#define VIPSTER_DEFINITIONS_H

#include "global.h"
#include <array>
#include <cmath>
#include <limits>
#include <algorithm>

namespace Vipster{

using Vec = std::array<float, 3>;
using Mat = std::array<Vec, 3>;

inline bool float_comp(float a, float b)
{
    return std::abs(a-b) <=
            std::numeric_limits<float>::epsilon() * 5 *
            std::max({1.f, std::abs(a), std::abs(b)});
}

inline bool operator==(const Vec &v1, const Vec &v2)
{
    return float_comp(v1[0], v2[0]) && float_comp(v1[1], v2[1]) && float_comp(v1[2], v2[2]);
}

inline bool operator!=(const Vec &v1, const Vec &v2)
{
    return !(v1==v2);
}

inline Vec& operator +=(Vec &v1, const Vec &v2)
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

inline Vec& operator +=(Vec &v, const float &f)
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

inline Vec& operator -=(Vec &v1, const Vec &v2)
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

inline Vec& operator -=(Vec &v, const float &f)
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

inline Vec& operator *=(Vec &v, const float &f)
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

inline Vec& operator /=(Vec &v, const float &f)
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

inline Mat& operator*=(Mat &m, const float &f)
{
    m[0]*=f;
    m[1]*=f;
    m[2]*=f;
    return m;
}

inline Mat operator*(Mat m, const float &f)
{
    return m*=f;
}

inline Mat operator*(const float &f, Mat m)
{
    return m*=f;
}

inline Mat& operator/=(Mat &m, const float&f)
{
    m[0]/=f;
    m[1]/=f;
    m[2]/=f;
    return m;
}

inline Mat operator/(Mat m, const float &f)
{
    return m/=f;
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

inline Mat& operator*=(Mat &lhs, const Mat &rhs)
{
    Mat tmp = lhs;
    Mat mt = Mat_trans(rhs);
    lhs[0][0] = Vec_dot(tmp[0],mt[0]);
    lhs[0][1] = Vec_dot(tmp[0],mt[1]);
    lhs[0][2] = Vec_dot(tmp[0],mt[2]);
    lhs[1][0] = Vec_dot(tmp[1],mt[0]);
    lhs[1][1] = Vec_dot(tmp[1],mt[1]);
    lhs[1][2] = Vec_dot(tmp[1],mt[2]);
    lhs[2][0] = Vec_dot(tmp[2],mt[0]);
    lhs[2][1] = Vec_dot(tmp[2],mt[1]);
    lhs[2][2] = Vec_dot(tmp[2],mt[2]);
    return lhs;
}

inline Mat operator*(Mat lhs, const Mat &rhs)
{
    return lhs*=rhs;
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
        throw Error("Mat_inv: singular matrix has no inverse!");
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

}

#endif // VIPSTER_DEFINITIONS_H
