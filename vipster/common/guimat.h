#ifndef GUIMAT_H
#define GUIMAT_H

#include <array>
#include "vec.h"

namespace Vipster {
namespace GUI {

typedef std::array<float,16> Mat;

}

GUI::Mat operator *=(GUI::Mat &a, const GUI::Mat &b);
GUI::Mat operator *(GUI::Mat a, const GUI::Mat &b);
void guiMatScale(GUI::Mat &m, float f);
void guiMatTranslate(GUI::Mat &m, float x, float y, float z);
void guiMatRot(GUI::Mat &m, float a, float x, float y, float z);
GUI::Mat guiMatMkOrtho(float left, float right, float bottom, float top, float near, float far);
GUI::Mat guiMatMkPerspective(float vertAng, float aspect, float near, float far);
GUI::Mat guiMatMkLookAt(Vec eye, Vec target, Vec up);

}

#endif // GUIMAT_H
