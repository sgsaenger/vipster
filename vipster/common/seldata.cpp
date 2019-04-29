#include "seldata.h"
#include "atom_model.h"

using namespace Vipster;

decltype(GUI::SelData::shader) GUI::SelData::shader;

GUI::SelData::SelData(const GUI::GlobalData& glob, const ColVec& col, Step::selection *sel)
    : Data{glob}, curSel{sel}
{
    color[0] = col[0]/255.f;
    color[1] = col[1]/255.f;
    color[2] = col[2]/255.f;
    color[3] = col[3]/255.f;
    update(curSel);
}

GUI::SelData::SelData(SelData&& dat)
//    : Data{dat.global, dat.updated, dat.initialized},
    : Data{std::move(dat)},
      sel_buffer{std::move(dat.sel_buffer)},
      cell_mat{dat.cell_mat},
      curSel{dat.curSel}
{
    std::swap(color, dat.color);
    std::swap(vao, dat.vao);
    std::swap(vbo, dat.vbo);
}

void GUI::SelData::initGL()
{
    if(!shader.initialized){
        shader.program = loadShader("/selection.vert", "/selection.frag");
        READATTRIB(shader, vertex)
        READATTRIB(shader, position)
        READATTRIB(shader, vert_scale)
        READATTRIB(shader, pbc_crit)
        READUNIFORM(shader, color)
        READUNIFORM(shader, offset)
        READUNIFORM(shader, pos_scale)
        READUNIFORM(shader, scale_fac)
        READUNIFORM(shader, mult)
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

    glVertexAttribIPointer(shader.pbc_crit, 3,
                          GL_SHORT,
                          sizeof(SelProp),
                          reinterpret_cast<const GLvoid*>(offsetof(SelProp, mult)));
    glVertexAttribDivisor(shader.pbc_crit, 1);
    glEnableVertexAttribArray(shader.pbc_crit);
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
        glUniform4fv(shader.color, 1, color);
        glUniform3i(shader.mult, 1, 1, 1);
        glDrawArraysInstanced(GL_TRIANGLES, 0,
                              atom_model_npoly,
                              static_cast<GLsizei>(sel_buffer.size()));
    }
}

void GUI::SelData::drawCell(const Vec &off, const PBCVec &mult)
{
    if(sel_buffer.size()){
        glBindVertexArray(vao);
        glUseProgram(shader.program);
        glUniform1f(shader.scale_fac, settings.atRadFac.val);
        glUniform3fv(shader.offset, 1, off.data());
        glUniformMatrix3fv(shader.pos_scale, 1, 0, cell_mat.data());
        glUniform4fv(shader.color, 1, color);
        glUniform3i(shader.mult, mult[0], mult[1], mult[2]);
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

void GUI::SelData::update(Step::selection* sel)
{
    curSel = sel;
    sel_buffer.clear();
    if(!curSel){
        return;
    }
    sel_buffer.reserve(curSel->getNat()); // too small, but better than nothing
    auto fmt = curSel->getFormatter(AtomFmt::Crystal, curSel->getFmt());
    if(settings.atRadVdW.val){
        auto it = curSel->begin();
        while(it != curSel->end()){
            for(const auto& off: it.getFilterPair().second){
                sel_buffer.push_back({it->coord + fmt(Vec{(float)off[0],(float)off[1],(float)off[2]}),
                                      it->pse->vdwr*1.3f,
                                      {(int16_t)off[0], (int16_t)off[1], (int16_t)off[2]},
                                     });
            }
            ++it;
        }
    }else{
        auto it = curSel->begin();
        while(it != curSel->end()){
            for(const auto& off: it.getFilterPair().second){
                sel_buffer.push_back({it->coord + fmt(Vec{(float)off[0],(float)off[1],(float)off[2]}),
                                      it->pse->covr*1.3f,
                                      {(int16_t)off[0], (int16_t)off[1], (int16_t)off[2]},
                                     });
            }
            ++it;
        }
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

void GUI::SelData::update(const ColVec &col)
{
    color[0] = col[0]/255.f;
    color[1] = col[1]/255.f;
    color[2] = col[2]/255.f;
    color[3] = col[3]/255.f;
}
