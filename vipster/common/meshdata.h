#ifndef MILLERDATA_H
#define MILLERDATA_H

#include "guidata.h"
#include "molecule.h"

namespace Vipster{
namespace GUI {

class MeshData: public Data{
    // CPU-Data:
    std::vector<Vec> vertices;
    Vec offset;
    std::array<GLuint, 4> color;
    StepProper* step;
    bool updated{true};
    // GPU-State/Data:
    GLuint vao{0}, vbo{0};
    // Shader:
    struct{
        GLuint program;
        GLuint vertex;
        GLint offset, color;
    }shader;
public:
    MeshData(const GlobalData& glob, std::vector<Vec>&& vertices,
             Vec offset, ColVec color, StepProper* step);
    MeshData(MeshData&& dat);
    ~MeshData() override;
    void drawMol() override;
    void drawCell(const std::array<uint8_t,3> &mult) override;
    void syncToGPU(void) override;
    void update(std::vector<Vec>&& vertices);
    void update(Vec offset);
    void update(const ColVec& color);
};

}
}

#endif // MILLERDATA_H
