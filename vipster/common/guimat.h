#ifndef GUIMAT_H
#define GUIMAT_H

#include <array>
#include "vec.h"

namespace Vipster {
namespace GUI {

typedef std::array<float,16> Mat_16f;

}

GUI::Mat_16f operator *=(GUI::Mat_16f &a, const GUI::Mat_16f &b);
GUI::Mat_16f operator *(GUI::Mat_16f a, const GUI::Mat_16f &b);
void guiMatScale(GUI::Mat_16f &m, float f);
void guiMatTranslate(GUI::Mat_16f &m, float x, float y, float z);
void guiMatRot(GUI::Mat_16f &m, float a, float x, float y, float z);
GUI::Mat_16f guiMatMkOrtho(float left, float right, float bottom, float top, float near, float far);
GUI::Mat_16f guiMatMkPerspective(float vertAng, float aspect, float near, float far);
GUI::Mat_16f guiMatMkLookAt(Vec eye, Vec target, Vec up);

}

#endif // GUIMAT_H
