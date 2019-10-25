#ifndef STEPDATA_H
#define STEPDATA_H

#include "guidata.h"
#include "molecule.h"

namespace Vipster{
namespace GUI {
// CPU-side buffers for render-data
struct AtomProp{ // 9 bytes + 3 bytes padding + 12 bytes directly from step
    float   rad;  // 4 bytes
    ColVec  col;  // 4 bytes
    uint8_t hide; // 1 byte
};
struct BondProp{ // 64 bytes
    float mat[9];        // 9*4 = 36 bytes
    Vec pos;             // 3*4 = 12 bytes
    int16_t mult[4];     // 4*2 = 8 bytes
    ColVec col_a, col_b; // 2*4 = 8 bytes
};

struct StepData: public Data{
    // CPU-Data:
    std::vector<AtomProp> atom_buffer{};
    std::vector<BondProp> bond_buffer{};
    std::array<Vec,8> cell_buffer{};
    std::array<float, 9>  cell_mat{};
    Step* curStep;
    // CPU-State:
    float atRadFac{0};
    // GPU-State:
    std::map<void*, GLuint[4]> vaos;
    // GPU-Data:
    bool vbo_initialized{false};
    GLuint vbos[4] = {0,0,0,0};
    GLuint &atom_prop_vbo{vbos[0]}, &atom_pos_vbo{vbos[1]};
    GLuint &bond_vbo{vbos[2]}, &cell_vbo{vbos[3]};
    // Shader:
    static struct{
        GLuint program;
        GLuint vertex, position, vert_scale, color, hide;
        GLint offset, pos_scale, scale_fac;
        bool initialized{false};
    } atom_shader;
    static struct{
        GLuint program;
        GLuint vertex, position, color1, color2, mMatrix, pbc_crit;
        GLint offset, pos_scale, pbc_cell, mult;
        bool initialized{false};
    } bond_shader;
    static struct{
        GLuint program;
        GLuint vertex;
        GLint offset;
        bool initialized{false};
    } cell_shader;
    static struct{
        GLuint program;
        GLuint vertex, position, vert_scale, hide;
        GLint offset, pos_scale, scale_fac, pbc_instance;
        bool initialized{false};
    } sel_shader;
    // Methods:
    StepData(const GlobalData& glob, Step* step=nullptr);
    ~StepData() override;
    StepData(StepData&& s);
    StepData(const StepData& s)=delete;
    StepData& operator=(const StepData& s)=delete;
    StepData& operator=(StepData&& s)=delete;
    void draw(const Vec &off, const PBCVec &mult, const Mat &cv,
              bool drawCell, void *context) override;
    void update(Step* step, bool useVdW, float atRadFac, float bondRad);
    void drawSel(Vec off, const PBCVec &mult, void *context);
private:
    void updateGL() override;
    void initGL(void *context) override;
    void initAtom(GLuint vao);
    void initBond(GLuint vao);
    void initCell(GLuint vao);
    void initSel(GLuint vao);
};
}
}

#endif // STEPDATA_H
