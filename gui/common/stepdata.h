#ifndef STEPDATA_H
#define STEPDATA_H

#include "guidata.h"
#include "vipster/molecule.h"

namespace Vipster{
namespace GUI {
// CPU-side buffers for render-data
struct AtomProp{ // 12 bytes + 9 bytes + 3 bytes padding
    float   pos[3];  // 3*4 = 12 bytes
    float   rad;  // 4 bytes
    ColVec  col;  // 4 bytes
    uint8_t hide; // 1 byte
};
struct BondProp{ // 64 bytes
    float mat[9];        // 9*4 = 36 bytes
    float pos[3];        // 3*4 = 12 bytes
    int16_t mult[4];     // 4*2 = 8 bytes
    ColVec col_a, col_b; // 2*4 = 8 bytes
};

struct StepData: public Data{
    // CPU-Data:
    std::vector<AtomProp> atom_buffer{};
    std::vector<BondProp> bond_buffer{};
    std::array<float, 24> cell_buffer{};
    std::array<float, 9>  cell_mat{};
    const Step* curStep;
    // CPU-State:
    float atRadFac{0};
    // GPU-State+Data (stored per instance and per context):
    struct ObjectContext{
        bool initialized{false};
        GLuint atom_vao{};
        GLuint bond_vao{};
        GLuint cell_vao{};
        GLuint sel_vao{};
        GLuint atom_vbo{};
        GLuint bond_vbo{};
        GLuint cell_vbo{};
    };
    std::map<void*, ObjectContext> object_map;
    // Shaders (stored per class and per context:
    struct ShaderContext{
        bool initialized{false};
        struct{
            GLuint program;
            GLuint vertex, position, vert_scale, color, hide;
            GLint offset, pos_scale, scale_fac;
        } atom_shader;
        struct{
            GLuint program;
            GLuint vertex, position, color1, color2, mMatrix, pbc_crit;
            GLint offset, pos_scale, pbc_cell, mult;
        } bond_shader;
        struct{
            GLuint program;
            GLuint vertex;
            GLint offset;
        } cell_shader;
        struct{
            GLuint program;
            GLuint vertex, position, vert_scale, hide;
            GLint offset, pos_scale, scale_fac, pbc_instance;
        } sel_shader;
    };
    static std::map<void*, ShaderContext> shader_map;
    // Methods:
    StepData(const Step* step=nullptr);
    ~StepData() override;
    StepData(StepData&& s);
    StepData(const StepData& s)=delete;
    StepData& operator=(const StepData& s)=delete;
    StepData& operator=(StepData&& s)=delete;
    void draw(const Vec &off, const PBCVec &mult, const Mat &cv,
              bool drawCell, void *context) override;
    void update(const Step* step, bool useVdW, float atRadFac, float bondRad);
    void drawSel(Vec off, const PBCVec &mult, void *context);
private:
    void updateGL(void *context) override;
    void initGL(void *context) override;
    void initAtom(GLuint vao);
    void initBond(GLuint vao);
    void initCell(GLuint vao);
    void initSel(GLuint vao);
    void initAtomVAO(GlobalContext& globals, ObjectContext& objects, ShaderContext& shaders);
    void initBondVAO(GlobalContext& globals, ObjectContext& objects, ShaderContext& shaders);
    void initCellVAO(GlobalContext& globals, ObjectContext& objects, ShaderContext& shaders);
    void initSelVAO(GlobalContext& globals, ObjectContext& objects, ShaderContext& shaders);
    void initAtomShader(GlobalContext& globals, ShaderContext& shaders);
    void initBondShader(GlobalContext& globals, ShaderContext& shaders);
    void initCellShader(GlobalContext& globals, ShaderContext& shaders);
    void initSelShader(GlobalContext& globals, ShaderContext& shaders);
};
}
}

#endif // STEPDATA_H
