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

class StepData: public Data{
    // CPU-Data:
    std::vector<AtomProp> atom_buffer{};
    std::vector<BondProp> bond_buffer{};
    std::array<Vec,8> cell_buffer{};
    std::array<float, 9>  cell_mat{};
    Step* curStep;
    // CPU-State:
    bool draw_bonds{false}, draw_cell{false};
    // GPU-State:
    GLuint vaos[4] = {0,0,0,0};
    GLuint &atom_vao{vaos[0]}, &bond_vao{vaos[1]};
    GLuint &cell_vao{vaos[2]}, &sel_vao{vaos[3]};
    // GPU-Data:
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
public:
    StepData(const GlobalData& glob, Step* step=nullptr);
    ~StepData() override;
    StepData(StepData&& s);
    StepData(const StepData& s)=delete;
    StepData& operator=(const StepData& s)=delete;
    StepData& operator=(StepData&& s)=delete;
    void drawMol(const Vec &off) override;
    void drawCell(const Vec &off, const PBCVec &mult) override;
    void updateGL() override;
    void initGL() override;
    void update(Step* step, bool b, bool c);
    void drawSel(const PBCVec &mult);
private:
    void initAtom();
    void initBond();
    void initCell();
    void initSel();
};
}
}

#endif // STEPDATA_H
