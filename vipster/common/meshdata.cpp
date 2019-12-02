#include "meshdata.h"

using namespace Vipster;

decltype (GUI::MeshData::mesh_shader) GUI::MeshData::mesh_shader;
decltype (GUI::MeshData::cell_shader) GUI::MeshData::cell_shader;

GUI::MeshData::MeshData(const GlobalData& glob, std::vector<Face>&& faces,
                        Vec offset, Mat cell, Texture texture)
    : Data{glob}, offset{offset},
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
//    : Data{dat.global, dat.updated, dat.initialized},
    : Data{std::move(dat)},
      offset{dat.offset},
      faces{std::move(dat.faces)},
      cell{dat.cell},
      cell_buffer{dat.cell_buffer},
      cell_gpu{dat.cell_gpu},
      texture{dat.texture}
{
    std::swap(vaos, dat.vaos);
    std::swap(vbos, dat.vbos);
    std::swap(tex, dat.tex);
}

void GUI::MeshData::initGL(void *context)
{
    auto &vao = vaos[context];
    glGenVertexArrays(2, vao);
    if(!vbo_initialized){
        glGenBuffers(2, vbos);
        glGenTextures(1, &tex);
        vbo_initialized = true;
    }

    initMesh(vao[0]);
    initCell(vao[1]);
    glBindVertexArray(0);
}

void GUI::MeshData::initMesh(GLuint vao)
{
    if(!mesh_shader.initialized){
        mesh_shader.program = loadShader("/mesh.vert", "/mesh.frag");
        READATTRIB(mesh_shader, vertex)
        READATTRIB(mesh_shader, normal)
        READATTRIB(mesh_shader, vert_UV)
        READUNIFORM(mesh_shader, pos_scale)
        READUNIFORM(mesh_shader, offset)
        READUNIFORM(mesh_shader, tex)
        mesh_shader.initialized = true;
    }

    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, mesh_vbo);
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

void GUI::MeshData::initCell(GLuint vao)
{
    if(!cell_shader.initialized){
        cell_shader.program = loadShader("/cell.vert", "/cell.frag");
        READATTRIB(cell_shader, vertex)
        READUNIFORM(cell_shader, offset)
        cell_shader.initialized = true;
    }

    glBindVertexArray(vao);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, global.cell_ibo);
    glBindBuffer(GL_ARRAY_BUFFER, cell_vbo);
    glVertexAttribPointer(cell_shader.vertex, 3, GL_FLOAT, GL_FALSE, 0, nullptr);
    glEnableVertexAttribArray(cell_shader.vertex);

}

GUI::MeshData::~MeshData()
{
    if(vbo_initialized){
        glDeleteBuffers(2, vbos);
        glDeleteTextures(1, &tex);
    }
    for(auto& vao: vaos){
        glDeleteVertexArrays(2, vao.second);
    }
}

void GUI::MeshData::update(std::vector<Face>&& face)
{
    faces = std::move(face);
    updated = true;
}

void GUI::MeshData::update(const Texture& tex)
{
    texture = tex;
    updated = true;
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
    updated = true;
}

void GUI::MeshData::updateGL()
{
    glBindBuffer(GL_ARRAY_BUFFER, mesh_vbo);
    glBufferData(GL_ARRAY_BUFFER, static_cast<GLsizeiptr>(sizeof(Face)*faces.size()),
                 static_cast<const void*>(faces.data()), GL_STREAM_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, cell_vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(cell_buffer),
                 static_cast<void*>(cell_buffer.data()), GL_STREAM_DRAW);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, tex);
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
    const auto& vao = vaos[context];
    glBindVertexArray(vao[0]);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, tex);
    glUseProgram(mesh_shader.program);
    glUniformMatrix3fv(mesh_shader.pos_scale, 1, 0, cell_gpu.data());
    glUniform1i(mesh_shader.tex, 0);
    for(int x=0;x<mult[0];++x){
        for(int y=0;y<mult[1];++y){
            for(int z=0;z<mult[2];++z){
                tmp2 = (tmp + x*cell[0] + y*cell[1] + z*cell[2]);
                glUniform3f(mesh_shader.offset,
                            static_cast<float>(tmp2[0]),
                            static_cast<float>(tmp2[1]),
                            static_cast<float>(tmp2[2]));
                glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(faces.size()));
            }
        }
    }
    if(drawCell){
        glBindVertexArray(vao[1]);
        glUseProgram(cell_shader.program);
        for(int x=0;x<mult[0];++x){
            for(int y=0;y<mult[1];++y){
                for(int z=0;z<mult[2];++z){
                    tmp2 = (tmp + x*cell[0] + y*cell[1] + z*cell[2]);
                    glUniform3f(mesh_shader.offset,
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
