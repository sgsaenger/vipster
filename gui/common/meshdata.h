#ifndef MILLERDATA_H
#define MILLERDATA_H

#include "guidata.h"
#include "vipster/global.h"

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
    struct ObjectContext{
        bool initialized{false};
        GLuint mesh_vao{};
        GLuint cell_vao{};
        GLuint mesh_vbo{};
        GLuint cell_vbo{};
        GLuint tex{};
    };
    std::map<void*, ObjectContext> object_map;
    // Shader:
    struct ShaderContext{
        struct{
            GLuint program;
            GLuint vertex, normal, vert_UV;
            GLint pos_scale, offset, tex;
        } mesh_shader;
        struct{
            GLuint program;
            GLuint vertex;
            GLint offset;
        } cell_shader;
        bool initialized{false};
    };
    static std::map<void*, ShaderContext> shader_map;
public:
    MeshData(std::vector<Face>&& faces, Vec offset,
             Mat cell, Texture texture);
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
    void updateGL(void *context) override;
    void initGL(void *context) override;
    void initMeshVAO(GlobalContext& globals, ObjectContext& objects, ShaderContext& shaders);
    void initCellVAO(GlobalContext& globals, ObjectContext& objects, ShaderContext& shaders);
    void initMeshShader(GlobalContext& globals, ShaderContext& shaders);
    void initCellShader(GlobalContext& globals, ShaderContext& shaders);
};

}
}

#endif // MILLERDATA_H
