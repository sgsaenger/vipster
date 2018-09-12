#include "meshdata.h"

using namespace Vipster;

GUI::MeshData::MeshData(const GlobalData& glob, std::vector<Vec>&& vertices,
                        Vec offset, Mat cell, ColVec color)
    : Data{glob}, vertices{std::move(vertices)},
      offset{offset}, cell{cell}, color{color}
{
    cell_gpu = {{cell[0][0], cell[1][0], cell[2][0],
                 cell[0][1], cell[1][1], cell[2][1],
                 cell[0][2], cell[1][2], cell[2][2]}};
}

GUI::MeshData::MeshData(MeshData&& dat)
    : Data{dat.global, dat.updated, dat.initialized},
      vertices{std::move(dat.vertices)},
      offset{dat.offset},
      cell{dat.cell},
      cell_gpu{dat.cell_gpu},
      color{dat.color}
{
    vao = dat.vao;
    dat.vao = 0;
    vbo = dat.vbo;
    dat.vbo = 0;
}

void GUI::MeshData::initGL()
{
    glGenVertexArrays(1, &vao);
    glGenBuffers(1, &vbo);

    shader.program = loadShader("/mesh.vert", "/mesh.frag");
    READATTRIB(shader, vertex);
    READUNIFORM(shader, pos_scale);
    READUNIFORM(shader, offset);
    READUNIFORM(shader, color);

    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glVertexAttribPointer(shader.vertex, 3, GL_FLOAT, GL_FALSE, 0, nullptr);
    glEnableVertexAttribArray(shader.vertex);
    glBindVertexArray(0);
}

GUI::MeshData::~MeshData()
{
    if(initialized){
        glDeleteBuffers(1, &vbo);
        glDeleteVertexArrays(1, &vao);
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

void GUI::MeshData::update(Vec off)
{
    offset = off;
}

void GUI::MeshData::update(Mat c)
{
    cell = c;
    cell_gpu = {{cell[0][0], cell[1][0], cell[2][0],
                 cell[0][1], cell[1][1], cell[2][1],
                 cell[0][2], cell[1][2], cell[2][2]}};
}

void GUI::MeshData::updateGL()
{
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, static_cast<GLsizeiptr>(sizeof(Vec)*vertices.size()),
                 static_cast<const void*>(vertices.data()), GL_STREAM_DRAW);
}

void GUI::MeshData::drawMol()
{
    glDisable(GL_CULL_FACE);
    Vec off = offset - (cell[0]+cell[1]+cell[2])/2;
    glBindVertexArray(vao);
    glUseProgram(shader.program);
    glUniformMatrix3fv(shader.pos_scale, 1, 0, cell_gpu.data());
    glUniform4f(shader.color, color[0]/255.f,
                              color[1]/255.f,
                              color[2]/255.f,
                              color[3]/255.f);
    glUniform3fv(shader.offset, 1, off.data());
    glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(vertices.size()));
    glEnable(GL_CULL_FACE);
}

void GUI::MeshData::drawCell(const std::array<uint8_t,3>& mult)
{
    glDisable(GL_CULL_FACE);
    Vec off;
    Vec tmp = offset - (cell[0]+cell[1]+cell[2])/2;
    tmp -= (mult[0]-1)*cell[0]/2.;
    tmp -= (mult[1]-1)*cell[1]/2.;
    tmp -= (mult[2]-1)*cell[2]/2.;
    glBindVertexArray(vao);
    glUseProgram(shader.program);
    glUniformMatrix3fv(shader.pos_scale, 1, 0, cell_gpu.data());
    glUniform4f(shader.color, color[0]/255.f,
                              color[1]/255.f,
                              color[2]/255.f,
                              color[3]/255.f);
    for(int x=0;x<mult[0];++x){
        for(int y=0;y<mult[1];++y){
            for(int z=0;z<mult[2];++z){
                off = (tmp + x*cell[0] + y*cell[1] + z*cell[2]);
                glUniform3fv(shader.offset, 1, off.data());
                glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(vertices.size()));
            }
        }
    }
    glEnable(GL_CULL_FACE);
}
