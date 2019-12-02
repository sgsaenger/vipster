#ifndef MILLERDATA_H
#define MILLERDATA_H

#include "guidata.h"
#include "global.h"

namespace Vipster{
namespace GUI {

class MeshData: public Data{
public:
    struct Face{
        float pos[3], norm[3];
        std::array<float, 2> uv;
    };
    struct Texture{
        std::vector<ColVec> data;
        int width, height;
    };
private:
    // CPU-Data:
public:
    Vec offset;
private:
    std::vector<Face> faces;
    Mat cell;
    std::array<float, 24> cell_buffer;
    std::array<float, 9> cell_gpu;
    Texture texture;
    // GPU-State/Data:
    std::map<void*, GLuint[2]> vaos;
    GLuint vbos[2] = {0,0};
    bool vbo_initialized{false};
    GLuint tex{0};
    GLuint &mesh_vbo{vbos[0]}, &cell_vbo{vbos[1]};
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
             Vec offset, Mat cell, Texture texture);
    MeshData(MeshData&& dat);
    MeshData(const MeshData& dat)=delete;
    MeshData& operator=(const MeshData& dat)=delete;
    MeshData& operator=(MeshData&& dat)=delete;
    ~MeshData() override;
    void draw(const Vec &off, const PBCVec &mult, const Mat &cv,
              bool drawCell, void *context) override;
    void update(std::vector<Face>&& faces);
    void update(const Texture& tex);
    void update(const Vipster::Mat &cell);
private:
    void updateGL(void) override;
    void initGL(void *context) override;
    void initMesh(GLuint vao);
    void initCell(GLuint vao);
};

}
}

#endif // MILLERDATA_H
