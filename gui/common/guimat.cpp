#include "guimat.h"
#include "vipster/global.h"

using namespace Vipster;

void Vipster::guiMatScale(GUI::Mat_16f &m, float f)
{
    for(size_t i=0;i<4;i++){
        m[i*4+0]*=f;
        m[i*4+1]*=f;
        m[i*4+2]*=f;
    }
}

void Vipster::guiMatTranslate(GUI::Mat_16f &m, float x, float y, float z)
{
    //assuming 0 0 0 1 in last row of m
    m[3]+=x;
    m[7]+=y;
    m[11]+=z;
}

void Vipster::guiMatRot(GUI::Mat_16f &m, float a, float x, float y, float z)
{
    if(float_comp(a,0.f)){
        return;
    }
    float tmp = a * static_cast<float>(deg2rad);
    float s = std::sin(tmp);
    float c = std::cos(tmp);
    if(float_comp(x,0.f)){
        if(float_comp(y,0.f)){
            if(!float_comp(z,0.f)){
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
        }else if(float_comp(z,0.f)){
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
    }else if(float_comp(y,0.f) && float_comp(z,0.f)){
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
        float ax[3] = {static_cast<float>(axis[0]),
                       static_cast<float>(axis[1]),
                       static_cast<float>(axis[2])};
        Vec axismc = axis*(1-c);
        float mc[3] = {static_cast<float>(axismc[0]),
                       static_cast<float>(axismc[1]),
                       static_cast<float>(axismc[2])};
        GUI::Mat_16f rotate{{c+mc[0]*ax[0], mc[1]*ax[0]-s*ax[2], mc[2]*ax[0]+s*ax[1], 0,
                             mc[0]*ax[1]+s*ax[2], c+mc[1]*ax[1], mc[2]*ax[1]-s*ax[0], 0,
                             mc[0]*ax[2]-s*ax[1], mc[1]*ax[2]+s*ax[0], c+mc[2]*ax[2], 0,
                             0,0,0,1}};
        m = rotate * m;
    }
}

GUI::Mat_16f Vipster::guiMatMkOrtho(float left, float right, float bottom, float top, float near, float far)
{
    return GUI::Mat_16f{{2/(right-left), 0, 0, (right+left)/(left-right),
                         0, 2/(top-bottom), 0, (top+bottom)/(bottom-top),
                         0, 0, (2/(near-far)), ((far+near)/(near-far)),
                         0, 0, 0, 1}};
}

GUI::Mat_16f Vipster::guiMatMkPerspective(float vertAng, float aspect, float near, float far)
{
    float rad = static_cast<float>(deg2rad) * vertAng/2;
    float sin = std::sin(rad);
    float cotan = std::cos(rad) / sin;
    float clip = far - near;
    return GUI::Mat_16f{{cotan/aspect, 0, 0, 0,
                     0, cotan, 0, 0,
                     0, 0, -(near + far)/clip, -(2*near*far)/clip,
                     0, 0, -1, 0}};
}

GUI::Mat_16f Vipster::guiMatMkLookAt(Vec eye, Vec target, Vec up)
{
    Vec dir = target - eye;
    dir /= Vec_length(dir);
    Vec r = Vec_cross(dir, up);
    r /= Vec_length(r);
    Vec u = Vec_cross(r, dir);
    return GUI::Mat_16f{{static_cast<float>(r[0]),
                         static_cast<float>(r[1]),
                         static_cast<float>(r[2]),
                         static_cast<float>(-Vec_dot(r, eye)),
                         static_cast<float>(u[0]),
                         static_cast<float>(u[1]),
                         static_cast<float>(u[2]),
                         static_cast<float>(-Vec_dot(u, eye)),
                         static_cast<float>(-dir[0]),
                         static_cast<float>(-dir[1]),
                         static_cast<float>(-dir[2]),
                         static_cast<float>(Vec_dot(dir, eye)),
                         0, 0, 0, 1}};
}

GUI::Mat_16f Vipster::operator *=(GUI::Mat_16f &a, const GUI::Mat_16f &b)
{
    a = GUI::Mat_16f{{a[0]*b[0]+a[1]*b[4]+a[2]*b[ 8]+a[3]*b[12],
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

GUI::Mat_16f Vipster::operator *(GUI::Mat_16f a, const GUI::Mat_16f &b)
{
    return a*=b;
}
