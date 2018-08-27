#include "stepdata.h"
#include "atom_model.h"
#include "bond_model.h"

#include <iostream>

using namespace Vipster;

GUI::StepData::StepData(GlobalData& glob, StepProper* step)
    : Data{glob},
      atom_vao{vaos[0]}, bond_vao{vaos[1]}, cell_vao{vaos[2]},
      atom_prop_vbo{vbos[0]}, atom_pos_vbo{vbos[1]},
      bond_vbo{vbos[2]}, cell_vbo{vbos[3]},
      curStep{step}
{
    glGenVertexArrays(3, vaos);
    glGenBuffers(4, vbos);
    GLuint loc;

// ATOMS
    glBindVertexArray(atom_vao);

    // sphere vertices
    glBindBuffer(GL_ARRAY_BUFFER, global.sphere_vbo);
    loc = static_cast<GLuint>(glGetAttribLocation(global.atom_program, "vertex_modelspace"));
    glVertexAttribPointer(loc, 3, GL_FLOAT, GL_FALSE, 0, nullptr);
    glEnableVertexAttribArray(loc);

    // atom positions
    glBindBuffer(GL_ARRAY_BUFFER, atom_pos_vbo);
    loc = static_cast<GLuint>(glGetAttribLocation(global.atom_program, "position_modelspace"));
    glVertexAttribPointer(loc, 3,
                          GL_FLOAT, GL_FALSE,
                          sizeof(Vec), nullptr);
    glVertexAttribDivisor(loc, 1);
    glEnableVertexAttribArray(loc);

    // atom properties
    glBindBuffer(GL_ARRAY_BUFFER, atom_prop_vbo);
    loc = static_cast<GLuint>(glGetAttribLocation(global.atom_program, "scale_modelspace"));
    glVertexAttribPointer(loc, 1,
                          GL_FLOAT, GL_FALSE,
                          sizeof(AtomProp),
                          reinterpret_cast<const GLvoid*>(offsetof(AtomProp, rad)));
    glVertexAttribDivisor(loc, 1);
    glEnableVertexAttribArray(loc);
    loc = static_cast<GLuint>(glGetAttribLocation(global.atom_program, "color_input"));
    glVertexAttribPointer(loc, 4,
                          GL_UNSIGNED_BYTE, GL_TRUE,
                          sizeof(AtomProp),
                          reinterpret_cast<const GLvoid*>(offsetof(AtomProp, col)));
    glVertexAttribDivisor(loc, 1);
    glEnableVertexAttribArray(loc);


// BONDS
    glBindVertexArray(bond_vao);

    // cylinder vertices
    glBindBuffer(GL_ARRAY_BUFFER, global.cylinder_vbo);
    loc = static_cast<GLuint>(glGetAttribLocation(global.bond_program, "vertex_modelspace"));
    glVertexAttribPointer(loc,3,GL_FLOAT,GL_FALSE,0,nullptr);
    glEnableVertexAttribArray(loc);

    // model matrix (rotation)
    glBindBuffer(GL_ARRAY_BUFFER, bond_vbo);
    loc = static_cast<GLuint>(glGetAttribLocation(global.bond_program, "mMatrix"));
    glVertexAttribPointer(loc, 3,
                          GL_FLOAT,GL_FALSE,
                          sizeof(BondProp), nullptr);
    glVertexAttribDivisor(loc, 1);
    glVertexAttribPointer(loc+1, 3,
                          GL_FLOAT,GL_FALSE,
                          sizeof(BondProp),
                          reinterpret_cast<const GLvoid*>(offsetof(BondProp, mat[3])));
    glVertexAttribDivisor(loc+1, 1);
    glVertexAttribPointer(loc+2, 3,
                          GL_FLOAT,GL_FALSE,
                          sizeof(BondProp),
                          reinterpret_cast<const GLvoid*>(offsetof(BondProp, mat[6])));
    glVertexAttribDivisor(loc+2, 1);
    glEnableVertexAttribArray(loc);
    glEnableVertexAttribArray(loc+1);
    glEnableVertexAttribArray(loc+2);

    // position
    loc = static_cast<GLuint>(glGetAttribLocation(global.bond_program, "position_modelspace"));
    glVertexAttribPointer(loc, 3,
                          GL_FLOAT,GL_FALSE,
                          sizeof(BondProp),
                          reinterpret_cast<const GLvoid*>(offsetof(BondProp, pos)));
    glVertexAttribDivisor(loc, 1);
    glEnableVertexAttribArray(loc);

    // pbc conditions
    loc = static_cast<GLuint>(glGetAttribLocation(global.bond_program, "pbc_crit"));
    glVertexAttribIPointer(loc, 4,
                          GL_UNSIGNED_SHORT,
                          sizeof(BondProp),
                          reinterpret_cast<const GLvoid*>(offsetof(BondProp, mult)));
    glVertexAttribDivisor(loc, 1);
    glEnableVertexAttribArray(loc);

    // Color 1
    loc = static_cast<GLuint>(glGetAttribLocation(global.bond_program, "s1Color"));
    glVertexAttribPointer(loc, 4,
                          GL_UNSIGNED_BYTE,GL_TRUE,
                          sizeof(BondProp),
                          reinterpret_cast<const GLvoid*>(offsetof(BondProp, col_a)));
    glVertexAttribDivisor(loc, 1);
    glEnableVertexAttribArray(loc);

    // Color 2
    loc = static_cast<GLuint>(glGetAttribLocation(global.bond_program, "s2Color"));
    glVertexAttribPointer(loc, 4,
                          GL_UNSIGNED_BYTE,GL_TRUE,
                          sizeof(BondProp),
                          reinterpret_cast<const GLvoid*>(offsetof(BondProp, col_b)));
    glVertexAttribDivisor(loc, 1);
    glEnableVertexAttribArray(loc);


// CELL
    glBindVertexArray(cell_vao);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, global.cell_ibo);
    glBindBuffer(GL_ARRAY_BUFFER, cell_vbo);
    loc = static_cast<GLuint>(glGetAttribLocation(global.cell_program, "vertex_modelspace"));
    glVertexAttribPointer(loc, 3, GL_FLOAT, GL_FALSE, 0, nullptr);
    glEnableVertexAttribArray(loc);
}

GUI::StepData::~StepData()
{
    glDeleteBuffers(4, vbos);
    glDeleteVertexArrays(3, vaos);
}

void GUI::StepData::drawMol()
{
    syncToGPU();
    Vec center = -curStep->getCenter(CdmFmt::Bohr);
// ATOMS
    glBindVertexArray(atom_vao);
    glUseProgram(global.atom_program);
    // TODO: move glGetUniformLocation to GlobalData and store locations
    auto offLocA = glGetUniformLocation(global.atom_program, "offset");
    auto facLoc = glGetUniformLocation(global.atom_program, "atom_fac");
    auto cellLocA = glGetUniformLocation(global.atom_program, "position_scale");
    auto toggleLoc = glGetUniformLocation(global.atom_program, "has_single_color");
    glUniform1f(facLoc, settings.atRadFac.val);
    glUniform3fv(offLocA, 1, center.data());
    glUniformMatrix3fv(cellLocA, 1, 0, cell_mat.data());
    glUniform1i(toggleLoc, 0);
    glDrawArraysInstanced(GL_TRIANGLES, 0,
                          atom_model_npoly,
                          static_cast<GLsizei>(atom_buffer.size()));
// BONDS
    if(settings.showBonds.val){
        glBindVertexArray(bond_vao);
        glUseProgram(global.bond_program);
        auto offLocB = glGetUniformLocation(global.bond_program, "offset");
        auto pbcLoc = glGetUniformLocation(global.bond_program, "pbc_cell");
        auto multLoc = glGetUniformLocation(global.bond_program, "mult");
        auto cellLocB = glGetUniformLocation(global.bond_program, "position_scale");
        glUniform3ui(multLoc, 1, 1, 1);
        glUniformMatrix3fv(cellLocB, 1, 0, cell_mat.data());
        glUniform3fv(offLocB, 1, center.data());
        glUniform3ui(pbcLoc, 0, 0, 0);
        glDrawArraysInstanced(GL_TRIANGLES, 0,
                              bond_model_npoly,
                              static_cast<GLsizei>(bond_buffer.size()));
    }
}

void GUI::StepData::drawCell(const std::array<uint8_t,3>& mult)
{
    syncToGPU();
    Vec off;
    Vec center = -curStep->getCenter(CdmFmt::Bohr);
    Mat cv = curStep->getCellVec() * curStep->getCellDim(CdmFmt::Bohr);
    center -= (mult[0]-1)*cv[0]/2.;
    center -= (mult[1]-1)*cv[1]/2.;
    center -= (mult[2]-1)*cv[2]/2.;
    // atoms
    glBindVertexArray(atom_vao);
    glUseProgram(global.atom_program);
    auto offLocA = glGetUniformLocation(global.atom_program, "offset");
    auto facLoc = glGetUniformLocation(global.atom_program, "atom_fac");
    auto cellLocA = glGetUniformLocation(global.atom_program, "position_scale");
    auto toggleLoc = glGetUniformLocation(global.atom_program, "has_single_color");
    glUniform1f(facLoc, settings.atRadFac.val);
    glUniform1f(toggleLoc, 0);
    glUniformMatrix3fv(cellLocA, 1, 0, cell_mat.data());
    for(int x=0;x<mult[0];++x){
        for(int y=0;y<mult[1];++y){
            for(int z=0;z<mult[2];++z){
                off = (center + x*cv[0] + y*cv[1] + z*cv[2]);
                glUniform3fv(offLocA, 1, off.data());
                glDrawArraysInstanced(GL_TRIANGLES, 0,
                                      atom_model_npoly,
                                      static_cast<GLsizei>(atom_buffer.size()));
            }
        }
    }
    // bonds
    if(settings.showBonds.val){
        glBindVertexArray(bond_vao);
        glUseProgram(global.bond_program);
        auto offLocB = glGetUniformLocation(global.bond_program, "offset");
        auto pbcLoc = glGetUniformLocation(global.bond_program, "pbc_cell");
        auto multLoc = glGetUniformLocation(global.bond_program, "mult");
        auto cellLocB = glGetUniformLocation(global.bond_program, "position_scale");
        glUniform3ui(multLoc, mult[0], mult[1], mult[2]);
        glUniformMatrix3fv(cellLocB, 1, 0, cell_mat.data());
        for(GLuint x=0;x<mult[0];++x){
            for(GLuint y=0;y<mult[1];++y){
                for(GLuint z=0;z<mult[2];++z){
                    off = (center + x*cv[0] + y*cv[1] + z*cv[2]);
                    glUniform3fv(offLocB, 1, off.data());
                    glUniform3ui(pbcLoc, x, y, z);
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
        glUseProgram(global.cell_program);
        auto offLocC = glGetUniformLocation(global.cell_program, "offset");
        for(int x=0;x<mult[0];++x){
            for(int y=0;y<mult[1];++y){
                for(int z=0;z<mult[2];++z){
                    off = (center + x*cv[0] + y*cv[1] + z*cv[2]);
                    glUniform3fv(offLocC, 1, off.data());
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
    syncToGPU();
#ifndef __EMSCRIPTEN__
    glDisable(GL_MULTISAMPLE);
#endif
    glClearColor(1,1,1,0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    Vec center = -curStep->getCenter(CdmFmt::Bohr);
    glBindVertexArray(atom_vao);
    glUseProgram(global.sel_program);
    auto offLocA = glGetUniformLocation(global.sel_program, "offset");
    auto facLoc = glGetUniformLocation(global.sel_program, "atom_fac");
    auto cellLocA = glGetUniformLocation(global.sel_program, "position_scale");
    glUniform1f(facLoc, settings.atRadFac.val);
    glUniformMatrix3fv(cellLocA, 1, 0, cell_mat.data());
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
                    glUniform3fv(offLocA, 1, off.data());
                    glDrawArraysInstanced(GL_TRIANGLES, 0,
                                          atom_model_npoly,
                                          static_cast<GLsizei>(atom_buffer.size()));
                }
            }
        }
    }else{
        glUniform3fv(offLocA, 1, center.data());
        glDrawArraysInstanced(GL_TRIANGLES, 0,
                              atom_model_npoly,
                              static_cast<GLsizei>(atom_buffer.size()));
    }
#ifndef __EMSCRIPTEN__
    glEnable(GL_MULTISAMPLE);
#endif
}

void GUI::StepData::syncToGPU()
{
    if (atoms_changed) {
        glBindBuffer(GL_ARRAY_BUFFER, atom_pos_vbo);
        auto nat = atom_buffer.size();
        if (nat != 0u) {
            glBufferData(GL_ARRAY_BUFFER, static_cast<GLsizeiptr>(nat*sizeof(Vec)),
                         static_cast<const void*>(curStep->getCoords().data()), GL_STREAM_DRAW);
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
        atoms_changed = false;
    }
    if (bonds_changed) {
        glBindBuffer(GL_ARRAY_BUFFER, bond_vbo);
        if (!bond_buffer.empty()) {
            glBufferData(GL_ARRAY_BUFFER, static_cast<GLsizeiptr>(bond_buffer.size()*sizeof(BondProp)),
                         static_cast<void*>(bond_buffer.data()), GL_STREAM_DRAW);
        } else {
            glBufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_STREAM_DRAW);
        }
        bonds_changed = false;
    }
    if (cell_changed) {
        glBindBuffer(GL_ARRAY_BUFFER, cell_vbo);
        glBufferData(GL_ARRAY_BUFFER, sizeof(cell_buffer),
                     static_cast<void*>(cell_buffer.data()), GL_STREAM_DRAW);
        cell_changed = false;
    }
}

void GUI::StepData::update(StepProper* step, bool draw_bonds)
{
    if (step!=nullptr) {
        curStep = step;
    }

// CELL
    Mat cv = curStep->getCellVec() * curStep->getCellDim(CdmFmt::Bohr);
    cell_buffer = {{ Vec{}, cv[0], cv[1], cv[2], cv[0]+cv[1], cv[0]+cv[2],
                     cv[1]+cv[2], cv[0]+cv[1]+cv[2] }};
    cell_changed = true;
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
    for (const PseEntry* at:curStep->getPseEntries()){
        atom_buffer.push_back({at->covr, at->col});
    }
    atoms_changed = true;

// BONDS
    if(draw_bonds){
        constexpr Vec x_axis{{1,0,0}};
        const auto& bonds = curStep->getBonds(settings.bondLvl.val);
        const auto& pse = curStep->getPseEntries();
        float c, s, ic;
        float rad = settings.bondRad.val;
        bond_buffer.clear();
        bond_buffer.reserve(bonds.size());
        const auto& at_coord = curStep->getCoords();
        auto fmt = curStep->getFmt();
        auto fmt_fun = curStep->getFormatter(fmt, AtomFmt::Bohr);
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
        bonds_changed = true;
        bonds_drawn = true;
    }else if(bonds_drawn){
        bond_buffer.clear();
        bonds_changed = true;
        bonds_drawn = false;
    }
}

bool GUI::StepData::hasCell() const noexcept
{
    return (curStep)?curStep->hasCell():false;
}
