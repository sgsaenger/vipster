#ifndef GLWRAPPER_H
#define GLWRAPPER_H

#include <string>
#include <GLES3/gl3.h>
#include <vector>
#include <array>
#include "molecule.h"

using namespace Vipster;

typedef std::array<float,16> glMat;
typedef std::array<float,4> glVec;
void glMatScale(glMat& m, float f);
void glMatRot(glMat& m, float a, float x, float y, float z);
//glMat glMatPerspective(float vAngle, float aspRatio, float near, float far);
glMat glMatOrtho(float left, float right, float bottom, float top, float near, float far);
glMat glMatLookAt(Vec eye, Vec target, Vec up);
glMat operator *=(glMat& a, const glMat& b);
glMat operator *(glMat a, const glMat& b);

GLuint loadShader(std::string header, std::string vertPath, std::string fragPath);

struct GLWrapper{
    GLWrapper();
    ~GLWrapper();
    // gpu-side data
    GLuint atom_program, bond_program, cell_program;
    GLuint atom_vao, bond_vao, cell_vao;
    GLuint atom_vbo, bond_vbo, cell_vbo;
    GLuint sphere_vbo, torus_vbo;
    // cpu-side data
    std::vector<std::array<float,8>> atom_buffer;
    std::vector<std::array<float,8>> bond_buffer;
    std::array<std::array<float,3>,8> cell_buffer;
    // cpu-side uniforms
    glMat vMat, pMat, rMat;
    bool vMatChanged, pMatChanged, rMatChanged;
    // other state
    int width, height;
};

#endif // GLWRAPPER_H
