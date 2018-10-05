#include "stepdata.h"
#include "atom_model.h"
#include "bond_model.h"

#include <iostream>

using namespace Vipster;

decltype(GUI::StepData::atom_shader) GUI::StepData::atom_shader;
decltype(GUI::StepData::bond_shader) GUI::StepData::bond_shader;
decltype(GUI::StepData::cell_shader) GUI::StepData::cell_shader;
decltype(GUI::StepData::sel_shader) GUI::StepData::sel_shader;

GUI::StepData::StepData(const GlobalData& glob, StepProper* step)
    : Data{glob},
      curStep{step}
{}

void GUI::StepData::initGL()
{
    glGenVertexArrays(4, vaos);
    glGenBuffers(4, vbos);

    initAtom();
    initBond();
    initCell();
    initSel();
    glBindVertexArray(0);
}

void GUI::StepData::initSel()
{
    if(!sel_shader.initialized){
        sel_shader.program = loadShader("/select.vert", "/select.frag");
        READATTRIB(sel_shader, vertex);
        READATTRIB(sel_shader, position);
        READATTRIB(sel_shader, vert_scale);
        READUNIFORM(sel_shader, pos_scale);
        READUNIFORM(sel_shader, scale_fac);
        READUNIFORM(sel_shader, offset);
        sel_shader.initialized = true;
    }

    glBindVertexArray(sel_vao);

    // sphere vertices
    glBindBuffer(GL_ARRAY_BUFFER, global.sphere_vbo);
    glVertexAttribPointer(sel_shader.vertex, 3, GL_FLOAT, GL_FALSE, 0, nullptr);
    glEnableVertexAttribArray(sel_shader.vertex);

    // atom positions
    glBindBuffer(GL_ARRAY_BUFFER, atom_pos_vbo);
    glVertexAttribPointer(sel_shader.position, 3,
                          GL_FLOAT, GL_FALSE,
                          sizeof(Vec), nullptr);
    glVertexAttribDivisor(sel_shader.position, 1);
    glEnableVertexAttribArray(sel_shader.position);

    // atom properties
    glBindBuffer(GL_ARRAY_BUFFER, atom_prop_vbo);
    glVertexAttribPointer(sel_shader.vert_scale, 1,
                          GL_FLOAT, GL_FALSE,
                          sizeof(AtomProp),
                          reinterpret_cast<const GLvoid*>(offsetof(AtomProp, rad)));
    glVertexAttribDivisor(sel_shader.vert_scale, 1);
    glEnableVertexAttribArray(sel_shader.vert_scale);
}

void GUI::StepData::initAtom()
{
    if(!atom_shader.initialized){
        atom_shader.program = loadShader("/atom.vert", "/atom.frag");
        READATTRIB(atom_shader, vertex);
        READATTRIB(atom_shader, position);
        READATTRIB(atom_shader, vert_scale);
        READATTRIB(atom_shader, color);
        READUNIFORM(atom_shader, offset);
        READUNIFORM(atom_shader, pos_scale);
        READUNIFORM(atom_shader, scale_fac);
        atom_shader.initialized = true;
    }

    glBindVertexArray(atom_vao);

    // sphere vertices
    glBindBuffer(GL_ARRAY_BUFFER, global.sphere_vbo);
    glVertexAttribPointer(atom_shader.vertex, 3, GL_FLOAT, GL_FALSE, 0, nullptr);
    glEnableVertexAttribArray(atom_shader.vertex);

    // atom positions
    glBindBuffer(GL_ARRAY_BUFFER, atom_pos_vbo);
    glVertexAttribPointer(atom_shader.position, 3,
                          GL_FLOAT, GL_FALSE,
                          sizeof(Vec), nullptr);
    glVertexAttribDivisor(atom_shader.position, 1);
    glEnableVertexAttribArray(atom_shader.position);

    // atom properties
    glBindBuffer(GL_ARRAY_BUFFER, atom_prop_vbo);
    glVertexAttribPointer(atom_shader.vert_scale, 1,
                          GL_FLOAT, GL_FALSE,
                          sizeof(AtomProp),
                          reinterpret_cast<const GLvoid*>(offsetof(AtomProp, rad)));
    glVertexAttribDivisor(atom_shader.vert_scale, 1);
    glEnableVertexAttribArray(atom_shader.vert_scale);
    glVertexAttribPointer(atom_shader.color, 4,
                          GL_UNSIGNED_BYTE, GL_TRUE,
                          sizeof(AtomProp),
                          reinterpret_cast<const GLvoid*>(offsetof(AtomProp, col)));
    glVertexAttribDivisor(atom_shader.color, 1);
    glEnableVertexAttribArray(atom_shader.color);
}

void GUI::StepData::initBond()
{
    if(!bond_shader.initialized){
        bond_shader.program = loadShader("/bond.vert", "/bond.frag");
        READATTRIB(bond_shader, vertex);
        READATTRIB(bond_shader, position);
        READATTRIB(bond_shader, color1);
        READATTRIB(bond_shader, color2);
        READATTRIB(bond_shader, mMatrix);
        READATTRIB(bond_shader, pbc_crit);
        READUNIFORM(bond_shader, offset);
        READUNIFORM(bond_shader, pos_scale);
        READUNIFORM(bond_shader, pbc_cell);
        READUNIFORM(bond_shader, mult);
        bond_shader.initialized = true;
    }

    glBindVertexArray(bond_vao);

    // cylinder vertices
    glBindBuffer(GL_ARRAY_BUFFER, global.cylinder_vbo);
    glVertexAttribPointer(bond_shader.vertex,3,GL_FLOAT,GL_FALSE,0,nullptr);
    glEnableVertexAttribArray(bond_shader.vertex);

    // model matrix (rotation)
    glBindBuffer(GL_ARRAY_BUFFER, bond_vbo);
    glVertexAttribPointer(bond_shader.mMatrix, 3,
                          GL_FLOAT,GL_FALSE,
                          sizeof(BondProp), nullptr);
    glVertexAttribDivisor(bond_shader.mMatrix, 1);
    glVertexAttribPointer(bond_shader.mMatrix+1, 3,
                          GL_FLOAT,GL_FALSE,
                          sizeof(BondProp),
                          reinterpret_cast<const GLvoid*>(offsetof(BondProp, mat[3])));
    glVertexAttribDivisor(bond_shader.mMatrix+1, 1);
    glVertexAttribPointer(bond_shader.mMatrix+2, 3,
                          GL_FLOAT,GL_FALSE,
                          sizeof(BondProp),
                          reinterpret_cast<const GLvoid*>(offsetof(BondProp, mat[6])));
    glVertexAttribDivisor(bond_shader.mMatrix+2, 1);
    glEnableVertexAttribArray(bond_shader.mMatrix);
    glEnableVertexAttribArray(bond_shader.mMatrix+1);
    glEnableVertexAttribArray(bond_shader.mMatrix+2);

    // position
    glVertexAttribPointer(bond_shader.position, 3,
                          GL_FLOAT,GL_FALSE,
                          sizeof(BondProp),
                          reinterpret_cast<const GLvoid*>(offsetof(BondProp, pos)));
    glVertexAttribDivisor(bond_shader.position, 1);
    glEnableVertexAttribArray(bond_shader.position);

    // pbc conditions
    glVertexAttribIPointer(bond_shader.pbc_crit, 4,
                          GL_UNSIGNED_SHORT,
                          sizeof(BondProp),
                          reinterpret_cast<const GLvoid*>(offsetof(BondProp, mult)));
    glVertexAttribDivisor(bond_shader.pbc_crit, 1);
    glEnableVertexAttribArray(bond_shader.pbc_crit);

    // Color 1
    glVertexAttribPointer(bond_shader.color1, 4,
                          GL_UNSIGNED_BYTE,GL_TRUE,
                          sizeof(BondProp),
                          reinterpret_cast<const GLvoid*>(offsetof(BondProp, col_a)));
    glVertexAttribDivisor(bond_shader.color1, 1);
    glEnableVertexAttribArray(bond_shader.color1);

    // Color 2
    glVertexAttribPointer(bond_shader.color2, 4,
                          GL_UNSIGNED_BYTE,GL_TRUE,
                          sizeof(BondProp),
                          reinterpret_cast<const GLvoid*>(offsetof(BondProp, col_b)));
    glVertexAttribDivisor(bond_shader.color2, 1);
    glEnableVertexAttribArray(bond_shader.color2);
}

void GUI::StepData::initCell()
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

GUI::StepData::~StepData()
{
    if(initialized){
        glDeleteBuffers(4, vbos);
        glDeleteVertexArrays(4, vaos);
    }
}

void GUI::StepData::drawMol(const Vec &off)
{
// ATOMS
    glBindVertexArray(atom_vao);
    glUseProgram(atom_shader.program);
    glUniform1f(atom_shader.scale_fac, settings.atRadFac.val);
    glUniform3fv(atom_shader.offset, 1, off.data());
    glUniformMatrix3fv(atom_shader.pos_scale, 1, 0, cell_mat.data());
    glDrawArraysInstanced(GL_TRIANGLES, 0,
                          atom_model_npoly,
                          static_cast<GLsizei>(atom_buffer.size()));
// BONDS
    if(settings.showBonds.val){
        glBindVertexArray(bond_vao);
        glUseProgram(bond_shader.program);
        glUniform3ui(bond_shader.mult, 1, 1, 1);
        glUniformMatrix3fv(bond_shader.pos_scale, 1, 0, cell_mat.data());
        glUniform3fv(bond_shader.offset, 1, off.data());
        glUniform3ui(bond_shader.pbc_cell, 0, 0, 0);
        glDrawArraysInstanced(GL_TRIANGLES, 0,
                              bond_model_npoly,
                              static_cast<GLsizei>(bond_buffer.size()));
    }
}

void GUI::StepData::drawCell(const Vec &off, const std::array<uint8_t,3>& mult)
{
    Vec tmp;
    Mat cv = curStep->getCellVec() * curStep->getCellDim(CdmFmt::Bohr);
    // atoms
    glBindVertexArray(atom_vao);
    glUseProgram(atom_shader.program);
    glUniform1f(atom_shader.scale_fac, settings.atRadFac.val);
    glUniformMatrix3fv(atom_shader.pos_scale, 1, 0, cell_mat.data());
    for(int x=0;x<mult[0];++x){
        for(int y=0;y<mult[1];++y){
            for(int z=0;z<mult[2];++z){
                tmp = (off + x*cv[0] + y*cv[1] + z*cv[2]);
                glUniform3fv(atom_shader.offset, 1, tmp.data());
                glDrawArraysInstanced(GL_TRIANGLES, 0,
                                      atom_model_npoly,
                                      static_cast<GLsizei>(atom_buffer.size()));
            }
        }
    }
    // bonds
    if(settings.showBonds.val){
        glBindVertexArray(bond_vao);
        glUseProgram(bond_shader.program);
        glUniform3ui(bond_shader.mult, mult[0], mult[1], mult[2]);
        glUniformMatrix3fv(bond_shader.pos_scale, 1, 0, cell_mat.data());
        for(GLuint x=0;x<mult[0];++x){
            for(GLuint y=0;y<mult[1];++y){
                for(GLuint z=0;z<mult[2];++z){
                    tmp = (off + x*cv[0] + y*cv[1] + z*cv[2]);
                    glUniform3fv(bond_shader.offset, 1, tmp.data());
                    glUniform3ui(bond_shader.pbc_cell, x, y, z);
                    glDrawArraysInstanced(GL_TRIANGLES, 0,
                                          bond_model_npoly,
                                          static_cast<GLsizei>(bond_buffer.size()));
                }
            }
        }
    }
    // cell
    if(settings.showCell.val){
        glBindVertexArray(cell_vao);
        glUseProgram(cell_shader.program);
        for(int x=0;x<mult[0];++x){
            for(int y=0;y<mult[1];++y){
                for(int z=0;z<mult[2];++z){
                    tmp = (off + x*cv[0] + y*cv[1] + z*cv[2]);
                    glUniform3fv(cell_shader.offset, 1, tmp.data());
                    glDrawElements(GL_LINES, 24, GL_UNSIGNED_SHORT, nullptr);
                }
            }
        }
    }
}

void GUI::StepData::drawSel(const std::array<uint8_t,3> &mult)
{
    // draw Atoms (as setup in atom_vao, shader locations must match)
    // with selection shader -> color by gl_InstanceID
    // to seperate Framebuffer
#ifndef __EMSCRIPTEN__
    if(settings.antialias.val){
        glDisable(GL_MULTISAMPLE);
    }
#endif
    glClearColor(1,1,1,0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    Vec center = -curStep->getCenter(CdmFmt::Bohr);
    glBindVertexArray(sel_vao);
    glUseProgram(sel_shader.program);
    glUniform1f(sel_shader.scale_fac, settings.atRadFac.val);
    glUniformMatrix3fv(sel_shader.pos_scale, 1, 0, cell_mat.data());
    if(curStep->hasCell()){
        Vec off;
        Mat cv = curStep->getCellVec() * curStep->getCellDim(CdmFmt::Bohr);
        center -= (mult[0]-1)*cv[0]/2.;
        center -= (mult[1]-1)*cv[1]/2.;
        center -= (mult[2]-1)*cv[2]/2.;
        for(int x=0;x<mult[0];++x){
            for(int y=0;y<mult[1];++y){
                for(int z=0;z<mult[2];++z){
                    off = (center + x*cv[0] + y*cv[1] + z*cv[2]);
                    glUniform3fv(sel_shader.offset, 1, off.data());
                    glDrawArraysInstanced(GL_TRIANGLES, 0,
                                          atom_model_npoly,
                                          static_cast<GLsizei>(atom_buffer.size()));
                }
            }
        }
    }else{
        glUniform3fv(atom_shader.offset, 1, center.data());
        glDrawArraysInstanced(GL_TRIANGLES, 0,
                              atom_model_npoly,
                              static_cast<GLsizei>(atom_buffer.size()));
    }
#ifndef __EMSCRIPTEN__
    if(settings.antialias.val){
        glEnable(GL_MULTISAMPLE);
    }
#endif
}

void GUI::StepData::updateGL()
{
    //TODO: separate data-handling somehow
    // ATOMS
    glBindBuffer(GL_ARRAY_BUFFER, atom_pos_vbo);
    auto nat = atom_buffer.size();
    if (nat != 0u) {
        glBufferData(GL_ARRAY_BUFFER, static_cast<GLsizeiptr>(nat*sizeof(Vec)),
                     static_cast<const void*>(
                     curStep->getAtoms().coordinates[static_cast<size_t>(curStep->getFmt())].data()
                     ), GL_STREAM_DRAW);
    } else {
        glBufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_STREAM_DRAW);
    }
    glBindBuffer(GL_ARRAY_BUFFER, atom_prop_vbo);
    if (nat != 0u) {
        glBufferData(GL_ARRAY_BUFFER, static_cast<GLsizeiptr>(nat*sizeof(AtomProp)),
                     static_cast<void*>(atom_buffer.data()), GL_STREAM_DRAW);
    } else {
        glBufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_STREAM_DRAW);
    }
    // BONDS
    glBindBuffer(GL_ARRAY_BUFFER, bond_vbo);
    if (!bond_buffer.empty()) {
        glBufferData(GL_ARRAY_BUFFER, static_cast<GLsizeiptr>(bond_buffer.size()*sizeof(BondProp)),
                     static_cast<void*>(bond_buffer.data()), GL_STREAM_DRAW);
    } else {
        glBufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_STREAM_DRAW);
    }
    // CELL
    glBindBuffer(GL_ARRAY_BUFFER, cell_vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(cell_buffer),
                 static_cast<void*>(cell_buffer.data()), GL_STREAM_DRAW);
}

void GUI::StepData::update(StepProper* step, bool draw_bonds)
{
    curStep = step;
    updated = true;

// CELL
    Mat cv = curStep->getCellVec() * curStep->getCellDim(CdmFmt::Bohr);
    cell_buffer = {{ Vec{}, cv[0], cv[1], cv[2], cv[0]+cv[1], cv[0]+cv[2],
                     cv[1]+cv[2], cv[0]+cv[1]+cv[2] }};
    Mat tmp_mat;
    if(curStep->getFmt() == AtomFmt::Crystal){
        tmp_mat = curStep->getCellVec();
    }else{
        tmp_mat = {{{{1,0,0}}, {{0,1,0}}, {{0,0,1}}}};
    }
    switch(curStep->getFmt()){
    case AtomFmt::Angstrom:
        tmp_mat *= Vipster::invbohr;
        break;
    case AtomFmt::Crystal:
    case AtomFmt::Alat:
        tmp_mat *= curStep->getCellDim(CdmFmt::Bohr);
        break;
    default:
        break;
    }
    cell_mat = {{tmp_mat[0][0], tmp_mat[1][0], tmp_mat[2][0],
                 tmp_mat[0][1], tmp_mat[1][1], tmp_mat[2][1],
                 tmp_mat[0][2], tmp_mat[1][2], tmp_mat[2][2]}};

// ATOMS
    atom_buffer.clear();
    atom_buffer.reserve(curStep->getNat());
    for (const auto& at: *curStep){
        atom_buffer.push_back({at.pse->covr, at.pse->col});
    }

// BONDS
    if(draw_bonds){
        constexpr Vec x_axis{{1,0,0}};
        const auto& bonds = curStep->getBonds();
        const auto& pse = curStep->getAtoms().pse;
        const auto& at_coord = curStep->getAtoms().coordinates[
                static_cast<size_t>(curStep->getFmt())];
        auto fmt = curStep->getFmt();
        auto fmt_fun = curStep->getFormatter(fmt, AtomFmt::Bohr);
        float c, s, ic;
        float rad = settings.bondRad.val;
        bond_buffer.clear();
        bond_buffer.reserve(bonds.size());
        switch(fmt){
        case AtomFmt::Crystal:
            cv = Mat{{{{1,0,0}}, {{0,1,0}}, {{0,0,1}}}};
            break;
        case AtomFmt::Alat:
            cv = curStep->getCellVec();
            break;
        case AtomFmt::Angstrom:
            cv = curStep->getCellVec() * curStep->getCellDim(CdmFmt::Angstrom);
            break;
        default:
            break;
        }
        Vec at_pos1, at_pos2, bond_pos, bond_axis, rot_axis;
        for(const Bond& bd:bonds){
            at_pos1 = at_coord[bd.at1];
            at_pos2 = at_coord[bd.at2];
            if (bd.xdiff>0)     { at_pos2 += bd.xdiff*cv[0]; }
            else if (bd.xdiff<0){ at_pos1 -= bd.xdiff*cv[0]; }
            if (bd.ydiff>0)     { at_pos2 += bd.ydiff*cv[1]; }
            else if (bd.ydiff<0){ at_pos1 -= bd.ydiff*cv[1]; }
            if (bd.zdiff>0)     { at_pos2 += bd.zdiff*cv[2]; }
            else if (bd.zdiff<0){ at_pos1 -= bd.zdiff*cv[2]; }
            bond_axis = at_pos1 - at_pos2;
            if(fmt == AtomFmt::Crystal){
                bond_axis = fmt_fun(bond_axis);
            }
            bond_pos = (at_pos1+at_pos2)/2;
            if(std::abs(bond_axis[1])<std::numeric_limits<float>::epsilon()&&
               std::abs(bond_axis[2])<std::numeric_limits<float>::epsilon()){
                c = std::copysign(1.f, bond_axis[0]);
                bond_buffer.push_back({
                    {bd.dist*c, 0., 0.,
                     0., rad, 0.,
                     0., 0., rad*c},
                    bond_pos,
                    {static_cast<uint16_t>(std::abs(bd.xdiff)),
                     static_cast<uint16_t>(std::abs(bd.ydiff)),
                     static_cast<uint16_t>(std::abs(bd.zdiff)),
                     static_cast<uint16_t>(!((bd.xdiff != 0)||(bd.ydiff != 0)||(bd.zdiff != 0)))},
                    pse[bd.at1]->col, pse[bd.at2]->col});
            }else{
                rot_axis = -Vec_cross(bond_axis, x_axis);
                rot_axis /= Vec_length(rot_axis);
                c = Vec_dot(bond_axis, x_axis)/Vec_length(bond_axis);
                ic = 1-c;
                s = -std::sqrt(1-c*c);
                bond_buffer.push_back({
                    //mat3 with rotation and scaling
                    {bd.dist*(ic*rot_axis[0]*rot_axis[0]+c),
                     bd.dist*(ic*rot_axis[0]*rot_axis[1]-s*rot_axis[2]),
                     bd.dist*(ic*rot_axis[0]*rot_axis[2]+s*rot_axis[1]),
                     rad*(ic*rot_axis[1]*rot_axis[0]+s*rot_axis[2]),
                     rad*(ic*rot_axis[1]*rot_axis[1]+c),
                     rad*(ic*rot_axis[1]*rot_axis[2]-s*rot_axis[0]),
                     rad*(ic*rot_axis[2]*rot_axis[0]-s*rot_axis[1]),
                     rad*(ic*rot_axis[2]*rot_axis[1]+s*rot_axis[0]),
                     rad*(ic*rot_axis[2]*rot_axis[2]+c)},
                    //vec3 with position in modelspace
                    bond_pos,
                    //faux uvec4 with integral pbc information
                    {static_cast<uint16_t>(std::abs(bd.xdiff)),
                     static_cast<uint16_t>(std::abs(bd.ydiff)),
                     static_cast<uint16_t>(std::abs(bd.zdiff)),
                    //padding that tells if non-pbc bond
                     static_cast<uint16_t>(!((bd.xdiff != 0)||(bd.ydiff != 0)||(bd.zdiff != 0)))},
                    //2*vec4 with colors
                    pse[bd.at1]->col, pse[bd.at2]->col});
            }
        }
        bonds_drawn = true;
    }else if(bonds_drawn){
        bond_buffer.clear();
        bonds_drawn = false;
    }
}
