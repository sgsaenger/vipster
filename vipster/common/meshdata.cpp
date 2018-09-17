#include "meshdata.h"

using namespace Vipster;

decltype (GUI::MeshData::mesh_shader) GUI::MeshData::mesh_shader;
decltype (GUI::MeshData::cell_shader) GUI::MeshData::cell_shader;

GUI::MeshData::MeshData(const GlobalData& glob, std::vector<Vec>&& vertices,
                        Vec offset, Mat cell, ColVec color)
    : Data{glob}, vertices{std::move(vertices)},
      offset{offset}, cell{cell}, color{color}
{
    cell_buffer = {{ Vec{}, cell[0], cell[1], cell[2], cell[0]+cell[1],
                   cell[0]+cell[2], cell[1]+cell[2], cell[0]+cell[1]+cell[2]}};
    cell_gpu = {{cell[0][0], cell[1][0], cell[2][0],
                 cell[0][1], cell[1][1], cell[2][1],
                 cell[0][2], cell[1][2], cell[2][2]}};
}

GUI::MeshData::MeshData(MeshData&& dat)
    : Data{dat.global, dat.updated, dat.initialized},
      vertices{std::move(dat.vertices)},
      offset{dat.offset},
      cell{dat.cell},
      cell_buffer{dat.cell_buffer},
      cell_gpu{dat.cell_gpu},
      color{dat.color}
{
    std::swap(mesh_vao, dat.mesh_vao);
    std::swap(mesh_vbo, dat.mesh_vbo);
    std::swap(cell_vao, dat.cell_vao);
    std::swap(cell_vbo, dat.cell_vbo);
}

void GUI::MeshData::initGL()
{
    glGenVertexArrays(2, vaos);
    glGenBuffers(2, vbos);

    initMesh();
    initCell();
    glBindVertexArray(0);
}

void GUI::MeshData::initMesh()
{
    if(!mesh_shader.initialized){
        mesh_shader.program = loadShader("/mesh.vert", "/mesh.frag");
        READATTRIB(mesh_shader, vertex);
        READUNIFORM(mesh_shader, pos_scale);
        READUNIFORM(mesh_shader, offset);
        READUNIFORM(mesh_shader, color);
        mesh_shader.initialized = true;
    }

    glBindVertexArray(mesh_vao);
    glBindBuffer(GL_ARRAY_BUFFER, mesh_vbo);
    glVertexAttribPointer(mesh_shader.vertex, 3, GL_FLOAT, GL_FALSE, 0, nullptr);
    glEnableVertexAttribArray(mesh_shader.vertex);
    glBindVertexArray(0);
}

void GUI::MeshData::initCell()
{
    if(!cell_shader.initialized){
        cell_shader.program = loadShader("/cell.vert", "/cell.frag");
        READATTRIB(cell_shader, vertex);
        READUNIFORM(cell_shader, offset);
        cell_shader.initialized = true;
    }

    glBindVertexArray(cell_vao);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, global.cell_ibo);
    glBindBuffer(GL_ARRAY_BUFFER, cell_vbo);
    glVertexAttribPointer(cell_shader.vertex, 3, GL_FLOAT, GL_FALSE, 0, nullptr);
    glEnableVertexAttribArray(cell_shader.vertex);

}

GUI::MeshData::~MeshData()
{
    if(initialized){
        glDeleteBuffers(2, vbos);
        glDeleteVertexArrays(2, vaos);
    }
}

void GUI::MeshData::update(std::vector<Vec>&& vert)
{
    vertices = std::move(vert);
    updated = true;
}

void GUI::MeshData::update(const ColVec& col)
{
    color = col;
}

void GUI::MeshData::update(Mat c)
{
    cell = c;
    cell_buffer = {{ Vec{}, cell[0], cell[1], cell[2], cell[0]+cell[1],
                   cell[0]+cell[2], cell[1]+cell[2], cell[0]+cell[1]+cell[2]}};
    cell_gpu = {{cell[0][0], cell[1][0], cell[2][0],
                 cell[0][1], cell[1][1], cell[2][1],
                 cell[0][2], cell[1][2], cell[2][2]}};
    updated = true;
}

void GUI::MeshData::updateGL()
{
    glBindBuffer(GL_ARRAY_BUFFER, mesh_vbo);
    glBufferData(GL_ARRAY_BUFFER, static_cast<GLsizeiptr>(sizeof(Vec)*vertices.size()),
                 static_cast<const void*>(vertices.data()), GL_STREAM_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, cell_vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(cell_buffer),
                 static_cast<void*>(cell_buffer.data()), GL_STREAM_DRAW);
}

void GUI::MeshData::drawMol(const Vec& off)
{
    glDisable(GL_CULL_FACE);
    Vec tmp = offset + off;
    glBindVertexArray(mesh_vao);
    glUseProgram(mesh_shader.program);
    glUniformMatrix3fv(mesh_shader.pos_scale, 1, 0, cell_gpu.data());
    glUniform4f(mesh_shader.color, color[0]/255.f,
                              color[1]/255.f,
                              color[2]/255.f,
                              color[3]/255.f);
    glUniform3fv(mesh_shader.offset, 1, tmp.data());
    glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(vertices.size()));
    glEnable(GL_CULL_FACE);
}

void GUI::MeshData::drawCell(const Vec& off, const std::array<uint8_t,3>& mult)
{
    glDisable(GL_CULL_FACE);
    Vec tmp = offset + off;
//    tmp -= (mult[0]-1)*cell[0]/2.;
//    tmp -= (mult[1]-1)*cell[1]/2.;
//    tmp -= (mult[2]-1)*cell[2]/2.;
    Vec tmp2;
    glBindVertexArray(mesh_vao);
    glUseProgram(mesh_shader.program);
    glUniformMatrix3fv(mesh_shader.pos_scale, 1, 0, cell_gpu.data());
    glUniform4f(mesh_shader.color, color[0]/255.f,
                              color[1]/255.f,
                              color[2]/255.f,
                              color[3]/255.f);
    for(int x=0;x<mult[0];++x){
        for(int y=0;y<mult[1];++y){
            for(int z=0;z<mult[2];++z){
                tmp2 = (tmp + x*cell[0] + y*cell[1] + z*cell[2]);
                glUniform3fv(mesh_shader.offset, 1, tmp2.data());
                glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(vertices.size()));
            }
        }
    }
    glBindVertexArray(cell_vao);
    glUseProgram(cell_shader.program);
    for(int x=0;x<mult[0];++x){
        for(int y=0;y<mult[1];++y){
            for(int z=0;z<mult[2];++z){
                tmp2 = (tmp + x*cell[0] + y*cell[1] + z*cell[2]);
                glUniform3fv(cell_shader.offset, 1, tmp2.data());
                glDrawElements(GL_LINES, 24, GL_UNSIGNED_SHORT, nullptr);
            }
        }
    }
    glEnable(GL_CULL_FACE);
}
