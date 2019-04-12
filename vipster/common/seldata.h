#ifndef SELDATA_H
#define SELDATA_H

#include "guidata.h"
#include "molecule.h"

namespace Vipster {
namespace GUI {
    struct SelProp{ // 22 bytes
        Vec pos;    // 3*4 = 12 bytes
        float rad;  // 4 bytes
        int16_t mult[3]; // 3*2 = 6 bytes
    };

    class SelData: public Data{
        // CPU-Data:
        std::vector<SelProp> sel_buffer{};
        std::array<float, 9>  cell_mat{};
        float color[4];
        Step::selection* curSel{nullptr};
        // GPU-State/Data:
        GLuint vao{0}, vbo{0};
        // Shader:
        static struct{
            GLuint program;
            GLuint vertex, position, vert_scale, pbc_crit;
            GLint offset, pos_scale, scale_fac, color, mult;
            bool initialized{false};
        } shader;
    public:
        SelData(const GlobalData& glob, const ColVec& color=settings.selCol.val,
                Step::selection* sel=nullptr);
        SelData(SelData&& dat);
        SelData& operator=(SelData&& dat)=delete;
        SelData(const SelData& dat)=delete;
        SelData& operator=(const SelData& dat)=delete;
        ~SelData() override;
        void drawMol(const Vec &off) override;
        void drawCell(const Vec &off, const PBCVec &mult) override;
        void updateGL() override;
        void initGL() override;
        void update(Step::selection* sel);
        void update(const ColVec& col);
    };
}
}

#endif // SELDATA_H
