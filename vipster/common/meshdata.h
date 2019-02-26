#ifndef MILLERDATA_H
#define MILLERDATA_H

#include "guidata.h"
#include "global.h"

namespace Vipster{
namespace GUI {

class MeshData: public Data{
public:
    struct Face{
        Vec pos, norm;
        std::array<float, 2> uv;
    };
    struct Texture{
        std::vector<ColVec> data;
        int width, height;
    };
private:
    // CPU-Data:
    std::vector<Face> faces;
    Vec offset;
    Mat cell;
    std::array<Vec, 8> cell_buffer;
    std::array<float, 9> cell_gpu;
    Texture texture;
    // GPU-State/Data:
    GLuint vaos[2] = {0,0};
    GLuint vbos[2] = {0,0};
    GLuint tex{0};
    GLuint &mesh_vao{vaos[0]}, &mesh_vbo{vbos[0]};
    GLuint &cell_vao{vaos[1]}, &cell_vbo{vbos[1]};
    // Shader:
    static struct{
        GLuint program;
        GLuint vertex, normal, vert_UV;
        GLint pos_scale, offset, tex;
        bool initialized{false};
    } mesh_shader;
    static struct{
        GLuint program;
        GLuint vertex;
        GLint offset;
        bool initialized{false};
    } cell_shader;
public:
    MeshData(const GlobalData& glob, std::vector<Face>&& faces,
             Vec offset, Vipster::Mat cell, Texture texture);
    MeshData(MeshData&& dat);
    MeshData(const MeshData& dat)=delete;
    MeshData& operator=(const MeshData& dat)=delete;
    MeshData& operator=(MeshData&& dat)=delete;
    ~MeshData() override;
    void drawMol(const Vec &off) override;
    void drawCell(const Vec &off, const PBCVec &mult) override;
    void updateGL(void) override;
    void initGL(void) override;
    void initMesh();
    void initCell();
    void update(std::vector<Face>&& faces);
    void update(const Texture& tex);
    void update(Vipster::Mat cell);
};

}
}

#endif // MILLERDATA_H
