#include "meshdata.h"

using namespace Vipster;

decltype (GUI::MeshData::shader_map) GUI::MeshData::shader_map;

GUI::MeshData::MeshData(std::vector<Face>&& faces,
                        Vec offset, Mat cell, Texture texture)
    : offset{offset},
      faces{std::move(faces)},
      cell{cell}, texture{texture}
{
    cell_buffer = {0,0,0,
                  static_cast<float>(cell[0][0]), static_cast<float>(cell[0][1]), static_cast<float>(cell[0][2]),
                  static_cast<float>(cell[1][0]), static_cast<float>(cell[1][1]), static_cast<float>(cell[1][2]),
                  static_cast<float>(cell[2][0]), static_cast<float>(cell[2][1]), static_cast<float>(cell[2][2]),
                  static_cast<float>(cell[0][0]+cell[1][0]), static_cast<float>(cell[0][1]+cell[1][1]), static_cast<float>(cell[0][2]+cell[1][2]),
                  static_cast<float>(cell[0][0]+cell[2][0]), static_cast<float>(cell[0][1]+cell[2][1]), static_cast<float>(cell[0][2]+cell[2][2]),
                  static_cast<float>(cell[1][0]+cell[2][0]), static_cast<float>(cell[1][1]+cell[2][1]), static_cast<float>(cell[1][2]+cell[2][2]),
                  static_cast<float>(cell[0][0]+cell[1][0]+cell[2][0]), static_cast<float>(cell[0][1]+cell[1][1]+cell[2][1]), static_cast<float>(cell[0][2]+cell[1][2]+cell[2][2]),
                  };
    cell_gpu = {{static_cast<float>(cell[0][0]),
                 static_cast<float>(cell[1][0]),
                 static_cast<float>(cell[2][0]),
                 static_cast<float>(cell[0][1]),
                 static_cast<float>(cell[1][1]),
                 static_cast<float>(cell[2][1]),
                 static_cast<float>(cell[0][2]),
                 static_cast<float>(cell[1][2]),
                 static_cast<float>(cell[2][2])}};
}

GUI::MeshData::MeshData(MeshData&& dat)
    : Data{std::move(dat)}
//      offset{dat.offset},
//      faces{std::move(dat.faces)},
//      cell{dat.cell},
//      cell_buffer{dat.cell_buffer},
//      cell_gpu{dat.cell_gpu},
//      texture{dat.texture}
{
    std::swap(offset, dat.offset);
    std::swap(faces, dat.faces);
    std::swap(cell, dat.cell);
    std::swap(cell_buffer, dat.cell_buffer);
    std::swap(cell_gpu, dat.cell_gpu);
    std::swap(texture, dat.texture);
    std::swap(object_map, dat.object_map);
}

void GUI::MeshData::initGL(void *context)
{
    auto &shaders = shader_map[context];
    if(!shaders.initialized){
        auto &globals = global_map[context];
        initMeshShader(globals, shaders);
        initCellShader(globals, shaders);
        shaders.initialized = true;
    }
    auto &objects = object_map[context];
    if(!objects.initialized){
        auto &globals = global_map[context];
        // VAOs
        glGenVertexArrays(1, &objects.mesh_vao);
        glGenVertexArrays(1, &objects.cell_vao);
        // VBOs
        glGenBuffers(1, &objects.mesh_vbo);
        glGenBuffers(1, &objects.cell_vbo);
        // texture
        glGenTextures(1, &objects.tex);
        // init
        initMeshVAO(globals, objects, shaders);
        initCellVAO(globals, objects, shaders);
        objects.initialized = true;
    }
    glBindVertexArray(0);
}

void GUI::MeshData::initMeshShader(GlobalContext& globals, ShaderContext &shaders)
{
    auto &mesh_shader = shaders.mesh_shader;
    mesh_shader.program = loadShader(globals, "/mesh.vert", "/mesh.frag");
    READATTRIB(mesh_shader, vertex)
    READATTRIB(mesh_shader, normal)
    READATTRIB(mesh_shader, vert_UV)
    READUNIFORM(mesh_shader, pos_scale)
    READUNIFORM(mesh_shader, offset)
    READUNIFORM(mesh_shader, tex)
}

void GUI::MeshData::initMeshVAO(GlobalContext&, ObjectContext &objects, ShaderContext &shaders)
{
    glBindVertexArray(objects.mesh_vao);
    glBindBuffer(GL_ARRAY_BUFFER, objects.mesh_vbo);

    auto &mesh_shader = shaders.mesh_shader;
    glVertexAttribPointer(mesh_shader.vertex, 3, GL_FLOAT, GL_FALSE, sizeof(Face),
                          reinterpret_cast<const GLvoid*>(offsetof(Face, pos)));
    glEnableVertexAttribArray(mesh_shader.vertex);
    glVertexAttribPointer(mesh_shader.normal, 3, GL_FLOAT, GL_FALSE, sizeof(Face),
                          reinterpret_cast<const GLvoid*>(offsetof(Face, norm)));
    glEnableVertexAttribArray(mesh_shader.normal);
    glVertexAttribPointer(mesh_shader.vert_UV, 2, GL_FLOAT, GL_TRUE, sizeof(Face),
                          reinterpret_cast<const GLvoid*>(offsetof(Face, uv)));
//    glVertexAttribIPointer(mesh_shader.vert_UV, 2, GL_UNSIGNED_BYTE, sizeof(Face),
//                          reinterpret_cast<const GLvoid*>(offsetof(Face, uv)));
    glEnableVertexAttribArray(mesh_shader.vert_UV);
}

void GUI::MeshData::initCellShader(GlobalContext& globals, ShaderContext &shaders)
{
    auto &cell_shader = shaders.cell_shader;
    cell_shader.program = loadShader(globals, "/cell.vert", "/cell.frag");
    READATTRIB(cell_shader, vertex)
    READUNIFORM(cell_shader, offset)
}

void GUI::MeshData::initCellVAO(GlobalContext& globals, ObjectContext &objects, ShaderContext &shaders)
{
    glBindVertexArray(objects.cell_vao);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, globals.cell_ibo);
    glBindBuffer(GL_ARRAY_BUFFER, objects.cell_vbo);
    glVertexAttribPointer(shaders.cell_shader.vertex, 3, GL_FLOAT, GL_FALSE, 0, nullptr);
    glEnableVertexAttribArray(shaders.cell_shader.vertex);

}

GUI::MeshData::~MeshData()
{
    for(auto &[context, objects]: object_map){
        if(!objects.initialized) continue;
        glDeleteVertexArrays(1, &objects.mesh_vao);
        glDeleteBuffers(1, &objects.mesh_vao);
        glDeleteVertexArrays(1, &objects.cell_vao);
        glDeleteBuffers(1, &objects.cell_vao);
        glDeleteTextures(1, &objects.tex);
    }
}

void GUI::MeshData::update(std::vector<Face>&& face)
{
    faces = std::move(face);
    for(auto &[context, state]: instance_map){
        state.synchronized = false;
    }
}

void GUI::MeshData::update(const Texture& tex)
{
    texture = tex;
    for(auto &[context, state]: instance_map){
        state.synchronized = false;
    }
}

void GUI::MeshData::update(const Mat &c)
{
    cell = c;
    cell_buffer = {0,0,0,
                  static_cast<float>(cell[0][0]), static_cast<float>(cell[0][1]), static_cast<float>(cell[0][2]),
                  static_cast<float>(cell[1][0]), static_cast<float>(cell[1][1]), static_cast<float>(cell[1][2]),
                  static_cast<float>(cell[2][0]), static_cast<float>(cell[2][1]), static_cast<float>(cell[2][2]),
                  static_cast<float>(cell[0][0]+cell[1][0]), static_cast<float>(cell[0][1]+cell[1][1]), static_cast<float>(cell[0][2]+cell[1][2]),
                  static_cast<float>(cell[0][0]+cell[2][0]), static_cast<float>(cell[0][1]+cell[2][1]), static_cast<float>(cell[0][2]+cell[2][2]),
                  static_cast<float>(cell[1][0]+cell[2][0]), static_cast<float>(cell[1][1]+cell[2][1]), static_cast<float>(cell[1][2]+cell[2][2]),
                  static_cast<float>(cell[0][0]+cell[1][0]+cell[2][0]), static_cast<float>(cell[0][1]+cell[1][1]+cell[2][1]), static_cast<float>(cell[0][2]+cell[1][2]+cell[2][2]),
                  };
    cell_gpu = {{static_cast<float>(cell[0][0]),
                 static_cast<float>(cell[1][0]),
                 static_cast<float>(cell[2][0]),
                 static_cast<float>(cell[0][1]),
                 static_cast<float>(cell[1][1]),
                 static_cast<float>(cell[2][1]),
                 static_cast<float>(cell[0][2]),
                 static_cast<float>(cell[1][2]),
                 static_cast<float>(cell[2][2])}};
    for(auto &[context, state]: instance_map){
        state.synchronized = false;
    }
}

void GUI::MeshData::updateGL(void *context)
{
    auto &objects = object_map[context];
    glBindBuffer(GL_ARRAY_BUFFER, objects.mesh_vbo);
    glBufferData(GL_ARRAY_BUFFER, static_cast<GLsizeiptr>(sizeof(Face)*faces.size()),
                 static_cast<const void*>(faces.data()), GL_STREAM_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, objects.cell_vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(cell_buffer),
                 static_cast<void*>(cell_buffer.data()), GL_STREAM_DRAW);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, objects.tex);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, texture.width, texture.height, 0, GL_RGBA,
                 GL_UNSIGNED_BYTE, static_cast<void*>(texture.data.data()));
    glBindTexture(GL_TEXTURE_2D, 0);
}

void GUI::MeshData::draw(const Vec &off, const PBCVec &mult,
                         const Mat &, bool drawCell, void *context)
{
    glDisable(GL_CULL_FACE);
    Vec tmp = offset + off;
    Vec tmp2;
    const auto& objects = object_map[context];
    const auto& shaders = shader_map[context];
    glBindVertexArray(objects.mesh_vao);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, objects.tex);
    glUseProgram(shaders.mesh_shader.program);
    glUniformMatrix3fv(shaders.mesh_shader.pos_scale, 1, 0, cell_gpu.data());
    glUniform1i(shaders.mesh_shader.tex, 0);
    for(int x=0;x<mult[0];++x){
        for(int y=0;y<mult[1];++y){
            for(int z=0;z<mult[2];++z){
                tmp2 = (tmp + x*cell[0] + y*cell[1] + z*cell[2]);
                glUniform3f(shaders.mesh_shader.offset,
                            static_cast<float>(tmp2[0]),
                            static_cast<float>(tmp2[1]),
                            static_cast<float>(tmp2[2]));
                glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(faces.size()));
            }
        }
    }
    if(drawCell){
        glBindVertexArray(objects.cell_vao);
        glUseProgram(shaders.cell_shader.program);
        for(int x=0;x<mult[0];++x){
            for(int y=0;y<mult[1];++y){
                for(int z=0;z<mult[2];++z){
                    tmp2 = (tmp + x*cell[0] + y*cell[1] + z*cell[2]);
                    glUniform3f(shaders.cell_shader.offset,
                                static_cast<float>(tmp2[0]),
                                static_cast<float>(tmp2[1]),
                                static_cast<float>(tmp2[2]));
                    glDrawElements(GL_LINES, 24, GL_UNSIGNED_SHORT, nullptr);
                }
            }
        }
    }
    glEnable(GL_CULL_FACE);
}
