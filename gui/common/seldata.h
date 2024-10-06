#ifndef SELDATA_H
#define SELDATA_H

#include "guidata.h"
#include "vipster/molecule.h"

namespace Vipster {
namespace GUI {
    struct SelProp{ // 22 bytes
        float pos[3];    // 3*4 = 12 bytes
        float rad;  // 4 bytes
        int16_t mult[3]; // 3*2 = 6 bytes
    };

    class SelData: public Data{
        // CPU-Data:
        std::vector<SelProp> sel_buffer{};
        std::array<float, 9>  cell_mat{};
        Step::const_selection const* curSel{nullptr};
        float atRadFac{};
        // GPU-State/Data:
        struct ObjectContext{
            bool initialized{false};
            GLuint vao{};
            GLuint vbo{};
        };
        std::map<void*, ObjectContext> object_map;
        // Shader:
        struct shader{
            GLuint program;
            GLuint vertex, position, vert_scale, pbc_crit;
            GLint offset, pos_scale, scale_fac, color, mult;
            bool initialized{false};
        };
        static std::map<void*, shader> shader_map;
        void updateGL(void *context) override;
        void initGL(void *context) override;
    public:
        ColVec color{};
        SelData(Step::const_selection const* sel=nullptr);
        SelData(SelData&& dat);
        SelData& operator=(SelData&& dat)=delete;
        SelData(const SelData& dat)=delete;
        SelData& operator=(const SelData& dat)=delete;
        ~SelData() override;
        void draw(const Vec &off, const PBCVec &mult, const Mat &cv,
                  bool drawCell, void *context) override;
        void update(Step::const_selection const* sel, bool useVdW, float atRadFac);
        void initShader(GlobalContext& globals, shader& shader);
        void initVAO(GlobalContext& globals, ObjectContext& objects, shader& shader);
    };
}
}

#endif // SELDATA_H
