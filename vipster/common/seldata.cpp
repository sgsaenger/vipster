#include "seldata.h"
#include "atom_model.h"

using namespace Vipster;

decltype(GUI::SelData::shader) GUI::SelData::shader;

GUI::SelData::SelData(const GUI::GlobalData& glob, StepSelection *sel)
    : Data{glob}, curSel{sel}
{}

void GUI::SelData::initGL()
{
    if(!shader.initialized){
        shader.program = loadShader("/selection.vert", "/selection.frag");
        READATTRIB(shader, vertex);
        READATTRIB(shader, position);
        READATTRIB(shader, vert_scale);
        READUNIFORM(shader, color);
        READUNIFORM(shader, offset);
        READUNIFORM(shader, pos_scale);
        READUNIFORM(shader, scale_fac);
        shader.initialized = true;
    }

    glGenVertexArrays(1, &vao);
    glGenBuffers(1, &vbo);
    glBindVertexArray(vao);

    glBindBuffer(GL_ARRAY_BUFFER, global.sphere_vbo);
    glVertexAttribPointer(shader.vertex, 3, GL_FLOAT, GL_FALSE, 0, nullptr);
    glEnableVertexAttribArray(shader.vertex);

    glBindBuffer(GL_ARRAY_BUFFER, vbo);

    glVertexAttribPointer(shader.position, 3,
                          GL_FLOAT, GL_FALSE,
                          sizeof(SelProp),
                          reinterpret_cast<const GLvoid*>(offsetof(SelProp, pos)));
    glVertexAttribDivisor(shader.position, 1);
    glEnableVertexAttribArray(shader.position);

    glVertexAttribPointer(shader.vert_scale, 1,
                          GL_FLOAT, GL_FALSE,
                          sizeof(SelProp),
                          reinterpret_cast<const GLvoid*>(offsetof(SelProp, rad)));
    glVertexAttribDivisor(shader.vert_scale, 1);
    glEnableVertexAttribArray(shader.vert_scale);
    glBindVertexArray(0);
}

GUI::SelData::~SelData()
{
    if(initialized){
        glDeleteBuffers(1, &vbo);
        glDeleteVertexArrays(1, &vao);
    }
}

void GUI::SelData::drawMol(const Vec &off)
{
    if(sel_buffer.size()){
        glBindVertexArray(vao);
        glUseProgram(shader.program);
        glUniform1f(shader.scale_fac, settings.atRadFac.val);
        glUniform3fv(shader.offset, 1, off.data());
        glUniformMatrix3fv(shader.pos_scale, 1, 0, cell_mat.data());
        glUniform4f(shader.color, settings.selCol.val[0]/255.f,
                                  settings.selCol.val[1]/255.f,
                                  settings.selCol.val[2]/255.f,
                                  settings.selCol.val[3]/255.f);
        glDrawArraysInstanced(GL_TRIANGLES, 0,
                              atom_model_npoly,
                              static_cast<GLsizei>(sel_buffer.size()));
    }
}

void GUI::SelData::drawCell(const Vec &off, const std::array<uint8_t, 3> &)
{
    if(sel_buffer.size()){
        glBindVertexArray(vao);
        glUseProgram(shader.program);
        glUniform1f(shader.scale_fac, settings.atRadFac.val);
        glUniform3fv(shader.offset, 1, off.data());
        glUniformMatrix3fv(shader.pos_scale, 1, 0, cell_mat.data());
        glUniform4f(shader.color, settings.selCol.val[0]/255.f,
                                  settings.selCol.val[1]/255.f,
                                  settings.selCol.val[2]/255.f,
                                  settings.selCol.val[3]/255.f);
        glDrawArraysInstanced(GL_TRIANGLES, 0,
                              atom_model_npoly,
                              static_cast<GLsizei>(sel_buffer.size()));
    }
}

void GUI::SelData::updateGL()
{
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    if(!sel_buffer.empty()){
        glBufferData(GL_ARRAY_BUFFER, static_cast<GLsizeiptr>(sel_buffer.size()*sizeof(SelProp)),
                     static_cast<const void*>(sel_buffer.data()), GL_STREAM_DRAW);
    }else{
        glBufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_STREAM_DRAW);
    }
}

void GUI::SelData::update(StepSelection* sel)
{
    curSel = sel;
    sel_buffer.clear();
    if(!curSel){
        return;
    }
    sel_buffer.reserve(curSel->getNat());
    for(const Atom& at:*curSel){
        sel_buffer.push_back({at.coord, at.pse->covr*1.3f});
    }
    updated = true;

    Mat tmp_mat;
    if(curSel->getFmt() == AtomFmt::Crystal){
        tmp_mat = curSel->getCellVec();
    }else{
        tmp_mat = {{{{1,0,0}}, {{0,1,0}}, {{0,0,1}}}};
    }
    switch(curSel->getFmt()){
    case AtomFmt::Angstrom:
        tmp_mat *= Vipster::invbohr;
        break;
    case AtomFmt::Crystal:
    case AtomFmt::Alat:
        tmp_mat *= curSel->getCellDim(CdmFmt::Bohr);
        break;
    default:
        break;
    }
    cell_mat = {{tmp_mat[0][0], tmp_mat[1][0], tmp_mat[2][0],
                 tmp_mat[0][1], tmp_mat[1][1], tmp_mat[2][1],
                 tmp_mat[0][2], tmp_mat[1][2], tmp_mat[2][2]}};
}
