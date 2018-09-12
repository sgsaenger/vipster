#ifndef MILLERDATA_H
#define MILLERDATA_H

#include "guidata.h"
#include "global.h"

namespace Vipster{
namespace GUI {

class MeshData: public Data{
    // CPU-Data:
    std::vector<Vec> vertices;
    Vec offset;
    Mat cell;
    std::array<float, 9> cell_gpu;
    std::array<uint8_t, 4> color;
    // GPU-State/Data:
    GLuint vao{0}, vbo{0};
    // Shader:
    struct{
        GLuint program;
        GLuint vertex;
        GLint pos_scale, offset, color;
    }shader;
public:
    MeshData(const GlobalData& glob, std::vector<Vec>&& vertices,
             Vec offset, Vipster::Mat cell, ColVec color);
    MeshData(MeshData&& dat);
    ~MeshData() override;
    void drawMol() override;
    void drawCell(const std::array<uint8_t,3> &mult) override;
    void updateGL(void) override;
    void initGL(void) override;
    void update(std::vector<Vec>&& vertices);
    void update(Vec offset);
    void update(const ColVec& color);
    void update(Vipster::Mat cell);
};

}
}

#endif // MILLERDATA_H
