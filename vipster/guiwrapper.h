#ifndef GLWRAPPER_H
#define GLWRAPPER_H

#include <string>
#ifdef __EMSCRIPTEN__
#include <GLES3/gl3.h>
#else
#include <QOpenGLFunctions_3_3_Core>
#endif
#include <vector>
#include <array>
#include "molecule.h"

using namespace Vipster;

typedef std::array<float,16> guiMat;

#ifdef __EMSCRIPTEN__
std::string readShader(std::string filePath);
struct GuiWrapper{
#else
#include <QString>
std::string readShader(QString filePath);
struct GuiWrapper: protected QOpenGLFunctions_3_3_Core{
#endif
    void loadShader(GLuint &program, std::string header, std::string vertShaderStr, std::string fragShaderStr);
    void initUBO(void);
    void initAtomVAO(void);
    void initBondVAO(void);
    void initCellVAO(void);
    void deleteGLObjects(void);
    // molecule-store
    std::vector<Vipster::Molecule> molecules;
    const Step* curStep{nullptr};
    // gpu-side data
    GLuint atom_program, bond_program, cell_program;
    GLuint atom_vao, bond_vao, cell_vao;
    GLuint atom_vbo, bond_vbo, cell_vbo;
    GLuint sphere_vbo, torus_vbo;
    GLuint cell_ibo;
    GLuint ubo;
    // cpu-side data
    std::vector<std::array<float,8>> atom_buffer;
    std::vector<std::array<float,24>> bond_buffer;
    std::array<Vec,8> cell_buffer;
    bool atoms_changed, bonds_changed, cell_changed;
    // cpu-side uniforms
    guiMat vMat, pMat, rMat;
    bool vMatChanged, pMatChanged, rMatChanged;
};

void guiMatScale(guiMat &m, float f);
void guiMatTranslate(guiMat &m, float x, float y, float z);
void guiMatRot(guiMat &m, float a, float x, float y, float z);
guiMat guiMatMkOrtho(float left, float right, float bottom, float top, float near, float far);
guiMat guiMatMkLookAt(Vec eye, Vec target, Vec up);
guiMat operator *=(guiMat &a, const guiMat &b);
guiMat operator *(guiMat a, const guiMat &b);

#endif // GLWRAPPER_H
