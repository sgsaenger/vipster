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

namespace Vipster {

typedef std::array<float,16> guiMat;

enum Change{atoms=0x1, cell=0x2, fmt=0x4, kpoints=0x8, param=0x10, config=0x20, selection=0x40};
constexpr auto stepChanged = Change::atoms | Change::cell | Change::fmt | Change::selection;
constexpr auto molChanged = Change::kpoints;

#ifdef __EMSCRIPTEN__
class GuiWrapper{
#else
class GuiWrapper: protected QOpenGLFunctions_3_3_Core{
#endif
    void loadShader(GLuint &program, const std::string &header, std::string vertShaderStr, std::string fragShaderStr);
public:
    void initShaders(const std::string& header, const std::string& folder);
    void deleteGLObjects(void);
    void draw(void);
    void drawCell(void);
    void drawMol(void);
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
    const StepProper* curStep{nullptr};
public:
    // cpu-side data
    std::array<uint8_t,3> mult{{1,1,1}};
private:
    // separate change-flag in step for coord and rest!
    struct atom_prop{ // 8 bytes + 12 bytes directly from step
        float rad;  // 4 bytes
        ColVec col; // 4 bytes
    };
    std::vector<atom_prop> atom_prop_buffer;
    struct bond_prop{ // 64 bytes
        float mat[9]; // 9*4 = 36 bytes
        Vec pos; // 3*4 = 12 bytes
        uint16_t mult[4];  // 4*2 = 6 bytes
        ColVec col_a, col_b; // 2*4 = 8 bytes
    };
    std::vector<bond_prop> bond_buffer{};
    std::array<Vec,8> cell_buffer{};
    std::array<float, 9>  cell_mat{};
    bool atoms_changed{false}, bonds_changed{false};
    bool cell_changed{false}, cell_mat_changed{false};
    bool bonds_drawn = false;
    // gpu-side data
    GLuint atom_program, bond_program, cell_program;
    GLuint atom_vao, bond_vao, cell_vao;
    GLuint atom_pos_vbo, atom_prop_vbo;
    GLuint bond_vbo, cell_vbo;
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

}

#endif // GLWRAPPER_H
