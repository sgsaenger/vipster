#ifndef STEPDATA_H
#define STEPDATA_H

#include "guidata.h"
#include "molecule.h"

namespace Vipster{
namespace GUI {
    // CPU-side buffers for render-data
    struct AtomProp{ // 8 bytes + 12 bytes directly from step
        float rad;  // 4 bytes
        ColVec col; // 4 bytes
    };
    struct BondProp{ // 64 bytes
        float mat[9]; // 9*4 = 36 bytes
        Vec pos; // 3*4 = 12 bytes
        uint16_t mult[4];  // 4*2 = 6 bytes
        ColVec col_a, col_b; // 2*4 = 8 bytes
    };

    class StepData: public Data{
        std::vector<AtomProp> atom_buffer{};
        std::vector<BondProp> bond_buffer{};
        std::array<Vec,8> cell_buffer{};
        std::array<float, 9>  cell_mat{};
        GLuint vaos[3];
        GLuint &atom_vao, &bond_vao, &cell_vao;
        GLuint vbos[4];
        GLuint &atom_prop_vbo, &atom_pos_vbo;
        GLuint &bond_vbo, &cell_vbo;
        StepProper* curStep;
        bool atoms_changed{false}, cell_changed{false};
        bool bonds_changed{false}, bonds_drawn{false};
    public:
        StepData(GlobalData& glob, StepProper* step=nullptr);
        ~StepData() override;
        void drawMol() override;
        void drawCell(const std::array<uint8_t,3> &mult) override;
        void syncToGPU() override;
        void update(StepProper* step, bool draw_bonds);
        void drawSel(const std::array<uint8_t,3> &mult);
        bool hasCell() const noexcept;
    };
}
}

#endif // STEPDATA_H
