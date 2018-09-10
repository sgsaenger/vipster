#include "meshdata.h"

using namespace Vipster;

GUI::MeshData::MeshData(const GlobalData& glob, std::vector<Vec>&& vertices,
                        Vec offset, ColVec color, StepProper* step)
    : Data{glob}, vertices{std::move(vertices)},
      offset{offset}, step{step}
{
    update(color);

    glGenVertexArrays(1, &vao);
    glGenBuffers(1, &vbo);

    shader.program = loadShader("/miller.vert", "/miller.frag");
    READATTRIB(shader, vertex);
    READUNIFORM(shader, offset);
    READUNIFORM(shader, color);

    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glVertexAttribPointer(shader.vertex, 3, GL_FLOAT, GL_FALSE, 0, nullptr);
    glEnableVertexAttribArray(shader.vertex);
    glBindVertexArray(0);
}

GUI::MeshData::MeshData(MeshData&& dat)
    : Data{dat.global},
      vertices{std::move(dat.vertices)},
      offset{dat.offset},
      color{dat.color},
      step{dat.step},
      updated{dat.updated}
{
    vao = dat.vao;
    dat.vao = 0;
    vbo = dat.vbo;
    dat.vbo = 0;
}

GUI::MeshData::~MeshData()
{
    glDeleteBuffers(1, &vbo);
    glDeleteVertexArrays(1, &vao);
}

void GUI::MeshData::update(std::vector<Vec>&& vert)
{
    vertices = std::move(vert);
    updated = true;
}

void GUI::MeshData::update(const ColVec& col)
{
    color[0] = col[0];
    color[1] = col[1];
    color[2] = col[2];
    color[3] = col[3];
}

void GUI::MeshData::update(Vec off)
{
    offset = off;
}

void GUI::MeshData::syncToGPU()
{
    if(updated){
        glBindBuffer(GL_ARRAY_BUFFER, vbo);
        glBufferData(GL_ARRAY_BUFFER, static_cast<GLsizeiptr>(sizeof(Vec)*vertices.size()),
                     static_cast<const void*>(vertices.data()), GL_STREAM_DRAW);
        updated = false;
    }
}

void GUI::MeshData::drawMol()
{
}

void GUI::MeshData::drawCell(const std::array<uint8_t,3>& mult)
{
    glDisable(GL_CULL_FACE);
    Vec off;
    Vec tmp = offset - step->getCenter(CdmFmt::Bohr);
    Mat cv = step->getCellVec() * step->getCellDim(CdmFmt::Bohr);
    tmp -= (mult[0]-1)*cv[0]/2.;
    tmp -= (mult[1]-1)*cv[1]/2.;
    tmp -= (mult[2]-1)*cv[2]/2.;
    glBindVertexArray(vao);
    glUseProgram(shader.program);
    glUniform4f(shader.color, color[0]/255.f,
                              color[1]/255.f,
                              color[2]/255.f,
                              color[3]/255.f);
    for(int x=0;x<mult[0];++x){
        for(int y=0;y<mult[1];++y){
            for(int z=0;z<mult[2];++z){
                off = (tmp + x*cv[0] + y*cv[1] + z*cv[2]);
                glUniform3fv(shader.offset, 1, off.data());
                glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(vertices.size()));
            }
        }
    }
    glEnable(GL_CULL_FACE);
}
