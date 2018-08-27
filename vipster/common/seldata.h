#ifndef SELDATA_H
#define SELDATA_H

#include "guidata.h"
#include "molecule.h"

namespace Vipster {
namespace GUI {
    struct SelProp{ // 16 bytes
        Vec pos;    // 3*4 = 12 bytes
        float rad;  // 4 bytes
    };

    class SelData: public Data{
        std::vector<SelProp> sel_buffer{};
        std::array<float, 9>  cell_mat{};
        GLuint vao, vbo;
        StepSelection* curSel;
        bool sel_changed{false};
    public:
        SelData(GlobalData& glob, StepSelection* sel=nullptr);
        ~SelData() override;
        void drawMol() override;
        void drawCell(const std::array<uint8_t,3> &mult) override;
        void syncToGPU() override;
        void update(StepSelection* sel);
    };
}
}

#endif // SELDATA_H
