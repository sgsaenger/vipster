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
class GuiWrapper{
#else
class GuiWrapper: protected QOpenGLFunctions_3_3_Core{
#endif
    void loadShader(GLuint &program, std::string header, std::string vertShaderStr, std::string fragShaderStr);
public:
    void initShaders(std::string header, std::string folder);
    void deleteGLObjects(void);
    void draw(void);
    // atom/bond/cell-data
    void initAtomVAO(void);
    void initBondVAO(void);
    void initCellVAO(void);
    void updateBuffers(const StepProper* step, bool draw_bonds=true);
    void updateVBOs(void);
    // view/projection matrices
    void initViewUBO(void);
    void updateViewUBO(void);
    void initViewMat(void);
    void resizeViewMat(int w, int h);
    void zoomViewMat(int i);
    void rotateViewMat(float x, float y, float z);
    void translateViewMat(float x, float y, float z);
    enum class alignDir{x,y,z,mx,my,mz};
    void alignViewMat(alignDir d);
    // molecule-store
    //TODO: is this needed here?
    std::vector<Vipster::Molecule> molecules;
    const StepProper* curStep{nullptr};
public:
    // cpu-side data
    std::array<int,3> mult{{1,1,1}};
private:
    std::vector<std::array<float,8>> atom_buffer;
    std::vector<std::array<float,24>> bond_buffer;
    std::array<Vec,8> cell_buffer;
    std::array<float, 9>  cell_mat;
    bool atoms_changed, bonds_changed, cell_changed, cell_mat_changed;
    bool bonds_drawn = false;
    // gpu-side data
    GLuint atom_program, bond_program, cell_program;
    GLuint atom_vao, bond_vao, cell_vao;
    GLuint atom_vbo, bond_vbo, cell_vbo;
    GLuint sphere_vbo, torus_vbo;
    GLuint cell_ibo;
    GLuint view_ubo;
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
