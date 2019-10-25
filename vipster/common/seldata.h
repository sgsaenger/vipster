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
        Step::selection* curSel{nullptr};
        // GPU-State/Data:
        float atRadFac{};
        std::map<void*, GLuint> vaos;
        GLuint vbo{0};
        bool vbo_initialized{false};
        // Shader:
        static struct{
            GLuint program;
            GLuint vertex, position, vert_scale, pbc_crit;
            GLint offset, pos_scale, scale_fac, color, mult;
            bool initialized{false};
        } shader;
        void updateGL() override;
        void initGL(void *context) override;
    public:
        ColVec color{};
        SelData(const GlobalData& glob, Step::selection* sel=nullptr);
        SelData(SelData&& dat);
        SelData& operator=(SelData&& dat)=delete;
        SelData(const SelData& dat)=delete;
        SelData& operator=(const SelData& dat)=delete;
        ~SelData() override;
        void draw(const Vec &off, const PBCVec &mult, const Mat &cv,
                  bool drawCell, void *context) override;
        void update(Step::selection* sel, bool useVdW, float atRadFac);
    };
}
}

#endif // SELDATA_H
