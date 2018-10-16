#include "guimat.h"
#include "global.h"

using namespace Vipster;

void Vipster::guiMatScale(GUI::Mat &m, float f)
{
    for(size_t i=0;i<4;i++){
        m[i*4+0]*=f;
        m[i*4+1]*=f;
        m[i*4+2]*=f;
    }
}

void Vipster::guiMatTranslate(GUI::Mat &m, float x, float y, float z)
{
    //assuming 0 0 0 1 in last row of m
    m[3]+=x;
    m[7]+=y;
    m[11]+=z;
}

void Vipster::guiMatRot(GUI::Mat &m, float a, float x, float y, float z)
{
    if(float_comp(a,0)){
        return;
    }
    float tmp = a * deg2rad;
    float s = std::sin(tmp);
    float c = std::cos(tmp);
    if(float_comp(x,0)){
        if(float_comp(y,0)){
            if(!float_comp(z,0)){
                // z-axis
                if (z<0){
                    s = -s;
                }
                m[0] = c * (tmp = m[0]) - s * m[4];
                m[4] = s * tmp + c * m[4];
                m[1] = c * (tmp = m[1]) - s * m[5];
                m[5] = s * tmp + c * m[5];
                m[2] = c * (tmp = m[2]) - s * m[6];
                m[6] = s * tmp + c * m[6];
                m[3] = c * (tmp = m[3]) - s * m[7];
                m[7] = s * tmp + c * m[7];
            }else{
                throw Error("guiMatRot: rotation axis-vector is zero");
            }
        }else if(float_comp(z,0)){
            // y-axis
            if (y<0) {
                s = -s;
            }
            m[0] = c * (tmp = m[0]) + s * m[8];
            m[8] = -s * tmp + c * m[8];
            m[1] = c * (tmp = m[1]) + s * m[9];
            m[9] = -s * tmp + c * m[9];
            m[2] = c * (tmp = m[2]) + s * m[10];
            m[10] = -s * tmp + c * m[10];
            m[3] = c * (tmp = m[3]) + s * m[11];
            m[11] = -s * tmp + c * m[11];
        }
    }else if(float_comp(y,0) && float_comp(z,0)){
        // x-axis
        if (x<0) {
            s = -s;
        }
        m[4] = c * (tmp = m[4]) - s * m[8];
        m[8] = s * tmp + c * m[8];
        m[5] = c * (tmp = m[5]) - s * m[9];
        m[9] = s * tmp + c * m[9];
        m[6] = c * (tmp = m[6]) - s * m[10];
        m[10] = s * tmp + c * m[10];
        m[7] = c * (tmp = m[7]) - s * m[11];
        m[11] = s * tmp + c * m[11];
    }else{
        // general rotation
        Vec axis{{x,y,z}};
        axis /= Vec_length(axis);
        Vec axismc = axis*(1-c);
        GUI::Mat rotate{{c+axismc[0]*axis[0], axismc[1]*axis[0]-s*axis[2], axismc[2]*axis[0]+s*axis[1], 0,
                      axismc[0]*axis[1]+s*axis[2], c+axismc[1]*axis[1], axismc[2]*axis[1]-s*axis[0], 0,
                      axismc[0]*axis[2]-s*axis[1], axismc[1]*axis[2]+s*axis[0], c+axismc[2]*axis[2], 0,
                      0,0,0,1}};
        m = rotate * m;
    }
}

GUI::Mat Vipster::guiMatMkOrtho(float left, float right, float bottom, float top, float near, float far)
{
    return GUI::Mat{{2/(right-left), 0, 0, (right+left)/(left-right),
                   0, 2/(top-bottom), 0, (top+bottom)/(bottom-top),
                   0, 0, (2/(near-far)), ((far+near)/(near-far)),
                   0, 0, 0, 1}};
}

GUI::Mat Vipster::guiMatMkPerspective(float vertAng, float aspect, float near, float far)
{
    float rad = deg2rad * vertAng/2;
    float sin = std::sin(rad);
    float cotan = std::cos(rad) / sin;
    float clip = far - near;
    return GUI::Mat{{cotan/aspect, 0, 0, 0,
                     0, cotan, 0, 0,
                     0, 0, -(near + far)/clip, -(2*near*far)/clip,
                     0, 0, -1, 0}};
}

GUI::Mat Vipster::guiMatMkLookAt(Vec eye, Vec target, Vec up)
{
    Vec dir = target - eye;
    dir /= Vec_length(dir);
    Vec r = Vec_cross(dir, up);
    r /= Vec_length(r);
    Vec u = Vec_cross(r, dir);
    return GUI::Mat{{r[0], r[1], r[2], -Vec_dot(r, eye),
                  u[0], u[1], u[2], -Vec_dot(u, eye),
                  -dir[0], -dir[1], -dir[2], Vec_dot(dir, eye),
                  0, 0, 0, 1}};
}

GUI::Mat Vipster::operator *=(GUI::Mat &a, const GUI::Mat &b)
{
    a = GUI::Mat{{a[0]*b[0]+a[1]*b[4]+a[2]*b[ 8]+a[3]*b[12],
               a[0]*b[1]+a[1]*b[5]+a[2]*b[ 9]+a[3]*b[13],
               a[0]*b[2]+a[1]*b[6]+a[2]*b[10]+a[3]*b[14],
               a[0]*b[3]+a[1]*b[7]+a[2]*b[11]+a[3]*b[15],
               a[4]*b[0]+a[5]*b[4]+a[6]*b[ 8]+a[7]*b[12],
               a[4]*b[1]+a[5]*b[5]+a[6]*b[ 9]+a[7]*b[13],
               a[4]*b[2]+a[5]*b[6]+a[6]*b[10]+a[7]*b[14],
               a[4]*b[3]+a[5]*b[7]+a[6]*b[11]+a[7]*b[15],
               a[8]*b[0]+a[9]*b[4]+a[10]*b[ 8]+a[11]*b[12],
               a[8]*b[1]+a[9]*b[5]+a[10]*b[ 9]+a[11]*b[13],
               a[8]*b[2]+a[9]*b[6]+a[10]*b[10]+a[11]*b[14],
               a[8]*b[3]+a[9]*b[7]+a[10]*b[11]+a[11]*b[15],
               a[12]*b[0]+a[13]*b[4]+a[14]*b[ 8]+a[15]*b[12],
               a[12]*b[1]+a[13]*b[5]+a[14]*b[ 9]+a[15]*b[13],
               a[12]*b[2]+a[13]*b[6]+a[14]*b[10]+a[15]*b[14],
               a[12]*b[3]+a[13]*b[7]+a[14]*b[11]+a[15]*b[15]}};
    return a;
}

GUI::Mat Vipster::operator *(GUI::Mat a, const GUI::Mat &b)
{
    return a*=b;
}
