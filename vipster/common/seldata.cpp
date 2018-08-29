#include "seldata.h"
#include "atom_model.h"

using namespace Vipster;

GUI::SelData::SelData(GUI::GlobalData& glob, StepSelection *sel)
    : Data{glob}, curSel{sel}
{
    glGenVertexArrays(1, &vao);
    glGenBuffers(1, &vbo);
    glBindVertexArray(vao);
    GLuint loc;
    glBindBuffer(GL_ARRAY_BUFFER, global.sphere_vbo);
    loc = static_cast<GLuint>(glGetAttribLocation(global.atom_program, "vertex_modelspace"));
    glVertexAttribPointer(loc, 3, GL_FLOAT, GL_FALSE, 0, nullptr);
    glEnableVertexAttribArray(loc);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);

    loc = static_cast<GLuint>(glGetAttribLocation(global.atom_program, "position_modelspace"));
    glVertexAttribPointer(loc, 3,
                          GL_FLOAT, GL_FALSE,
                          sizeof(SelProp),
                          reinterpret_cast<const GLvoid*>(offsetof(SelProp, pos)));
    glVertexAttribDivisor(loc, 1);
    glEnableVertexAttribArray(loc);

    loc = static_cast<GLuint>(glGetAttribLocation(global.atom_program, "scale_modelspace"));
    glVertexAttribPointer(loc, 1,
                          GL_FLOAT, GL_FALSE,
                          sizeof(SelProp),
                          reinterpret_cast<const GLvoid*>(offsetof(SelProp, rad)));
    glVertexAttribDivisor(loc, 1);
    glEnableVertexAttribArray(loc);
}

GUI::SelData::~SelData()
{
    glDeleteBuffers(1, &vbo);
    glDeleteVertexArrays(1, &vao);
}

void GUI::SelData::drawMol()
{
    if(sel_buffer.size()){
        Vec center = -curSel->getCenter(CdmFmt::Bohr);
        auto offLocA = glGetUniformLocation(global.atom_program, "offset");
        auto facLoc = glGetUniformLocation(global.atom_program, "atom_fac");
        auto cellLocA = glGetUniformLocation(global.atom_program, "position_scale");
        auto toggleLoc = glGetUniformLocation(global.atom_program, "has_single_color");
        auto colorLoc = glGetUniformLocation(global.atom_program, "single_color");
        glBindVertexArray(vao);
        glUseProgram(global.atom_program);
        glUniform1f(facLoc, settings.atRadFac.val);
        glUniform3fv(offLocA, 1, center.data());
        glUniformMatrix3fv(cellLocA, 1, 0, cell_mat.data());
        glUniform1i(toggleLoc, 1);
        glUniform4f(colorLoc, settings.selCol.val[0]/255.f,
                              settings.selCol.val[1]/255.f,
                              settings.selCol.val[2]/255.f,
                              settings.selCol.val[3]/255.f);
        glDrawArraysInstanced(GL_TRIANGLES, 0,
                              atom_model_npoly,
                              static_cast<GLsizei>(sel_buffer.size()));
    }
}

void GUI::SelData::drawCell(const std::array<uint8_t, 3> &mult)
{
    syncToGPU();
    if(sel_buffer.size()){
        Vec center = -curSel->getCenter(CdmFmt::Bohr);
        Mat cv = curSel->getCellVec() * curSel->getCellDim(CdmFmt::Bohr);
        center -= (mult[0]-1)*cv[0]/2.;
        center -= (mult[1]-1)*cv[1]/2.;
        center -= (mult[2]-1)*cv[2]/2.;
        auto offLocA = glGetUniformLocation(global.atom_program, "offset");
        auto facLoc = glGetUniformLocation(global.atom_program, "atom_fac");
        auto cellLocA = glGetUniformLocation(global.atom_program, "position_scale");
        auto toggleLoc = glGetUniformLocation(global.atom_program, "has_single_color");
        auto colorLoc = glGetUniformLocation(global.atom_program, "single_color");
        glBindVertexArray(vao);
        glUseProgram(global.atom_program);
        glUniform1f(facLoc, settings.atRadFac.val);
        glUniform3fv(offLocA, 1, center.data());
        glUniformMatrix3fv(cellLocA, 1, 0, cell_mat.data());
        glUniform1i(toggleLoc, 1);
        glUniform4f(colorLoc, settings.selCol.val[0]/255.f,
                              settings.selCol.val[1]/255.f,
                              settings.selCol.val[2]/255.f,
                              settings.selCol.val[3]/255.f);
        glDrawArraysInstanced(GL_TRIANGLES, 0,
                              atom_model_npoly,
                              static_cast<GLsizei>(sel_buffer.size()));
    }
}

void GUI::SelData::syncToGPU()
{
    if (sel_changed) {
        glBindBuffer(GL_ARRAY_BUFFER, vbo);
        if(!sel_buffer.empty()){
            glBufferData(GL_ARRAY_BUFFER, static_cast<GLsizeiptr>(sel_buffer.size()*sizeof(SelProp)),
                         static_cast<const void*>(sel_buffer.data()), GL_STREAM_DRAW);
        }else{
            glBufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_STREAM_DRAW);
        }
        sel_changed = false;
    }
}

void GUI::SelData::update(StepSelection* sel)
{
    if(sel != nullptr){
        curSel = sel;
    }
    sel_buffer.clear();
    sel_buffer.reserve(curSel->getNat());
    for(const Atom& at:*curSel){
        sel_buffer.push_back({at.coord, at.pse->covr*1.3f});
    }
    sel_changed = true;

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
