#include "stepdata.h"
#include "atom_model.h"
#include "bond_model.h"

using namespace Vipster;

decltype(GUI::StepData::atom_shader) GUI::StepData::atom_shader;
decltype(GUI::StepData::bond_shader) GUI::StepData::bond_shader;
decltype(GUI::StepData::cell_shader) GUI::StepData::cell_shader;
decltype(GUI::StepData::sel_shader) GUI::StepData::sel_shader;

// TODO: update itself if applicable if scaling can be moved to shader
GUI::StepData::StepData(const GlobalData& glob, Step* step)
    : Data{glob},
      curStep{step}
{}

GUI::StepData::StepData(StepData&& dat)
    : Data{std::move(dat)},
      curStep{dat.curStep}
{
    std::swap(vaos, dat.vaos);
    std::swap(vbos, dat.vbos);
}

GUI::StepData::~StepData()
{
    if(vbo_initialized){
        glDeleteBuffers(3, vbos);
    }
    for(auto &vao: vaos){
        glDeleteVertexArrays(4, vao.second);
    }
}

void GUI::StepData::initGL(void *context)
{
    auto &vao = vaos[context];
    glGenVertexArrays(4, vao);
    if(!vbo_initialized){
        glGenBuffers(3, vbos);
        vbo_initialized = true;
    }

    initAtom(vao[0]);
    initBond(vao[1]);
    initCell(vao[2]);
    initSel(vao[3]);
    glBindVertexArray(0);
}

void GUI::StepData::initSel(GLuint vao)
{
    if(!sel_shader.initialized){
        sel_shader.program = loadShader("/select.vert", "/select.frag");
        READATTRIB(sel_shader, vertex)
        READATTRIB(sel_shader, position)
        READATTRIB(sel_shader, vert_scale)
        READATTRIB(sel_shader, hide)
        READUNIFORM(sel_shader, pos_scale)
        READUNIFORM(sel_shader, scale_fac)
        READUNIFORM(sel_shader, offset)
        READUNIFORM(sel_shader, pbc_instance)
        sel_shader.initialized = true;
    }

    glBindVertexArray(vao);

    // sphere vertices
    glBindBuffer(GL_ARRAY_BUFFER, global.sphere_vbo);
    glVertexAttribPointer(sel_shader.vertex, 3, GL_FLOAT, GL_FALSE, 0, nullptr);
    glEnableVertexAttribArray(sel_shader.vertex);

    // atom positions
    glBindBuffer(GL_ARRAY_BUFFER, atom_vbo);
    glVertexAttribPointer(sel_shader.position, 3,
                          GL_FLOAT, GL_FALSE,
                          sizeof(AtomProp), nullptr);
    glVertexAttribDivisor(sel_shader.position, 1);
    glEnableVertexAttribArray(sel_shader.position);

    // atom properties
    glBindBuffer(GL_ARRAY_BUFFER, atom_vbo);
    glVertexAttribPointer(sel_shader.vert_scale, 1,
                          GL_FLOAT, GL_FALSE,
                          sizeof(AtomProp),
                          reinterpret_cast<const GLvoid*>(offsetof(AtomProp, rad)));
    glVertexAttribDivisor(sel_shader.vert_scale, 1);
    glEnableVertexAttribArray(sel_shader.vert_scale);
    glVertexAttribIPointer(sel_shader.hide, 1,
                           GL_UNSIGNED_BYTE,
                           sizeof(AtomProp),
                           reinterpret_cast<const GLvoid*>(offsetof(AtomProp, hide)));
    glVertexAttribDivisor(sel_shader.hide, 1);
    glEnableVertexAttribArray(sel_shader.hide);
}

void GUI::StepData::initAtom(GLuint vao)
{
    if(!atom_shader.initialized){
        atom_shader.program = loadShader("/atom.vert", "/atom.frag");
        READATTRIB(atom_shader, vertex)
        READATTRIB(atom_shader, position)
        READATTRIB(atom_shader, vert_scale)
        READATTRIB(atom_shader, color)
        READATTRIB(atom_shader, hide)
        READUNIFORM(atom_shader, offset)
        READUNIFORM(atom_shader, pos_scale)
        READUNIFORM(atom_shader, scale_fac)
        atom_shader.initialized = true;
    }

    glBindVertexArray(vao);

    // sphere vertices
    glBindBuffer(GL_ARRAY_BUFFER, global.sphere_vbo);
    glVertexAttribPointer(atom_shader.vertex, 3, GL_FLOAT, GL_FALSE, 0, nullptr);
    glEnableVertexAttribArray(atom_shader.vertex);

    glBindBuffer(GL_ARRAY_BUFFER, atom_vbo);
    // atom positions
    glVertexAttribPointer(atom_shader.position, 3,
                          GL_FLOAT, GL_FALSE,
                          sizeof(AtomProp),
                          reinterpret_cast<const GLvoid*>(offsetof(AtomProp, pos)));
    glVertexAttribDivisor(atom_shader.position, 1);
    glEnableVertexAttribArray(atom_shader.position);

    // atom properties
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
    glVertexAttribIPointer(atom_shader.hide, 1,
                           GL_UNSIGNED_BYTE,
                           sizeof(AtomProp),
                           reinterpret_cast<const GLvoid*>(offsetof(AtomProp, hide)));
    glVertexAttribDivisor(atom_shader.hide, 1);
    glEnableVertexAttribArray(atom_shader.hide);
}

void GUI::StepData::initBond(GLuint vao)
{
    if(!bond_shader.initialized){
        bond_shader.program = loadShader("/bond.vert", "/bond.frag");
        READATTRIB(bond_shader, vertex)
        READATTRIB(bond_shader, position)
        READATTRIB(bond_shader, color1)
        READATTRIB(bond_shader, color2)
        READATTRIB(bond_shader, mMatrix)
        READATTRIB(bond_shader, pbc_crit)
        READUNIFORM(bond_shader, offset)
        READUNIFORM(bond_shader, pos_scale)
        READUNIFORM(bond_shader, pbc_cell)
        READUNIFORM(bond_shader, mult)
        bond_shader.initialized = true;
    }

    glBindVertexArray(vao);

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
                           GL_SHORT,
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

void GUI::StepData::initCell(GLuint vao)
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

void GUI::StepData::draw(const Vec &off, const PBCVec &mult,
                         const Mat &cv, bool drawCell, void *context)
{
    Vec tmp;
    const auto& vao = vaos[context];
    // atoms
    glBindVertexArray(vao[0]);
    glUseProgram(atom_shader.program);
    glUniform1f(atom_shader.scale_fac, atRadFac);
    glUniformMatrix3fv(atom_shader.pos_scale, 1, 0, cell_mat.data());
    for(int x=0;x<mult[0];++x){
        for(int y=0;y<mult[1];++y){
            for(int z=0;z<mult[2];++z){
                tmp = (off + x*cv[0] + y*cv[1] + z*cv[2]);
                glUniform3f(atom_shader.offset,
                            static_cast<float>(tmp[0]),
                            static_cast<float>(tmp[1]),
                            static_cast<float>(tmp[2]));
                glDrawArraysInstanced(GL_TRIANGLES, 0,
                                      atom_model_npoly,
                                      static_cast<GLsizei>(atom_buffer.size()));
            }
        }
    }
    // bonds
    glBindVertexArray(vao[1]);
    glUseProgram(bond_shader.program);
    glUniform3i(bond_shader.mult, mult[0], mult[1], mult[2]);
    glUniformMatrix3fv(bond_shader.pos_scale, 1, 0, cell_mat.data());
    for(GLint x=0;x<mult[0];++x){
        for(GLint y=0;y<mult[1];++y){
            for(GLint z=0;z<mult[2];++z){
                tmp = (off + x*cv[0] + y*cv[1] + z*cv[2]);
                glUniform3f(bond_shader.offset,
                            static_cast<float>(tmp[0]),
                            static_cast<float>(tmp[1]),
                            static_cast<float>(tmp[2]));
                glUniform3i(bond_shader.pbc_cell, x, y, z);
                glDrawArraysInstanced(GL_TRIANGLES, 0,
                                      bond_model_npoly,
                                      static_cast<GLsizei>(bond_buffer.size()));
            }
        }
    }
    // cell
    if(drawCell && curStep->hasCell()){
        glBindVertexArray(vao[2]);
        glUseProgram(cell_shader.program);
        for(int x=0;x<mult[0];++x){
            for(int y=0;y<mult[1];++y){
                for(int z=0;z<mult[2];++z){
                    tmp = (off + x*cv[0] + y*cv[1] + z*cv[2]);
                    glUniform3f(cell_shader.offset,
                                static_cast<float>(tmp[0]),
                                static_cast<float>(tmp[1]),
                                static_cast<float>(tmp[2]));
                    glDrawElements(GL_LINES, 24, GL_UNSIGNED_SHORT, nullptr);
                }
            }
        }
    }
}

void GUI::StepData::drawSel(Vec off, const PBCVec &mult, void *context)
{
    // draw Atoms (as setup in atom_vao, shader locations must match)
    // with selection shader -> color by gl_InstanceID
    // to seperate Framebuffer
    glEnable(GL_CULL_FACE);
    glEnable(GL_DEPTH_TEST);
    glClearColor(1,1,0,0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glBindVertexArray(vaos[context][3]);
    glUseProgram(sel_shader.program);
    glUniform1f(sel_shader.scale_fac, atRadFac);
    glUniformMatrix3fv(sel_shader.pos_scale, 1, 0, cell_mat.data());
    if(curStep->hasCell()){
        Vec tmp;
        Mat cv = curStep->getCellVec() * curStep->getCellDim(AtomFmt::Bohr);
        off -= (mult[0]-1)*cv[0]/2.;
        off -= (mult[1]-1)*cv[1]/2.;
        off -= (mult[2]-1)*cv[2]/2.;
        for(unsigned int x=0;x<mult[0];++x){
            for(unsigned int y=0;y<mult[1];++y){
                for(unsigned int z=0;z<mult[2];++z){
                    tmp = (off + x*cv[0] + y*cv[1] + z*cv[2]);
                    glUniform3f(sel_shader.offset,
                                static_cast<float>(tmp[0]),
                                static_cast<float>(tmp[1]),
                                static_cast<float>(tmp[2]));
                    glUniform1ui(sel_shader.pbc_instance,
                                 1 + x + y*mult[0] + z*mult[0]*mult[1]);
                    glDrawArraysInstanced(GL_TRIANGLES, 0,
                                          atom_model_npoly,
                                          static_cast<GLsizei>(atom_buffer.size()));
                }
            }
        }
    }else{
        glUniform3f(sel_shader.offset,
                    static_cast<float>(off[0]),
                    static_cast<float>(off[1]),
                    static_cast<float>(off[2]));
        glUniform1ui(sel_shader.pbc_instance, 1);
        glDrawArraysInstanced(GL_TRIANGLES, 0,
                              atom_model_npoly,
                              static_cast<GLsizei>(atom_buffer.size()));
    }
}

void GUI::StepData::updateGL()
{
    //TODO: separate data-handling somehow
    // ATOMS
    auto nat = atom_buffer.size();
    glBindBuffer(GL_ARRAY_BUFFER, atom_vbo);
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

// TODO: radii should probably be factored in in shader
void GUI::StepData::update(Step* step,
                           bool useVdW, float atRadFac,
                           float bondRad)
{
    curStep = step;
    updated = true;
    this->atRadFac = atRadFac;

// CELL
    Mat cv = curStep->getCellVec() * curStep->getCellDim(AtomFmt::Bohr);
    cell_buffer = {0,0,0,
                  static_cast<float>(cv[0][0]), static_cast<float>(cv[0][1]), static_cast<float>(cv[0][2]),
                  static_cast<float>(cv[1][0]), static_cast<float>(cv[1][1]), static_cast<float>(cv[1][2]),
                  static_cast<float>(cv[2][0]), static_cast<float>(cv[2][1]), static_cast<float>(cv[2][2]),
                  static_cast<float>(cv[0][0]+cv[1][0]), static_cast<float>(cv[0][1]+cv[1][1]), static_cast<float>(cv[0][2]+cv[1][2]),
                  static_cast<float>(cv[0][0]+cv[2][0]), static_cast<float>(cv[0][1]+cv[2][1]), static_cast<float>(cv[0][2]+cv[2][2]),
                  static_cast<float>(cv[1][0]+cv[2][0]), static_cast<float>(cv[1][1]+cv[2][1]), static_cast<float>(cv[1][2]+cv[2][2]),
                  static_cast<float>(cv[0][0]+cv[1][0]+cv[2][0]), static_cast<float>(cv[0][1]+cv[1][1]+cv[2][1]), static_cast<float>(cv[0][2]+cv[1][2]+cv[2][2]),
                  };
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
        tmp_mat *= curStep->getCellDim(AtomFmt::Bohr);
        break;
    default:
        break;
    }
    cell_mat = {{static_cast<float>(tmp_mat[0][0]),
                 static_cast<float>(tmp_mat[1][0]),
                 static_cast<float>(tmp_mat[2][0]),
                 static_cast<float>(tmp_mat[0][1]),
                 static_cast<float>(tmp_mat[1][1]),
                 static_cast<float>(tmp_mat[2][1]),
                 static_cast<float>(tmp_mat[0][2]),
                 static_cast<float>(tmp_mat[1][2]),
                 static_cast<float>(tmp_mat[2][2])}};

// ATOMS
    atom_buffer.clear();
    atom_buffer.reserve(curStep->getNat());
    if(useVdW){
        for (const auto& at: *curStep){
            atom_buffer.push_back({{static_cast<float>(at.coord[0]),
                                    static_cast<float>(at.coord[1]),
                                    static_cast<float>(at.coord[2])},
                                   static_cast<float>(at.type->vdwr), at.type->col,
                                   static_cast<uint8_t>(at.properties->flags[AtomFlag::Hidden])});
        }
    }else{
        for (const auto& at: *curStep){
            atom_buffer.push_back({{static_cast<float>(at.coord[0]),
                                    static_cast<float>(at.coord[1]),
                                    static_cast<float>(at.coord[2])},
                                   static_cast<float>(at.type->covr), at.type->col,
                                   static_cast<uint8_t>(at.properties->flags[AtomFlag::Hidden])});
        }
    }

// BONDS
    constexpr Vec x_axis{{1,0,0}};
    const auto& bonds = curStep->getBonds();
    const auto& elements = curStep->getAtoms().elements;
    const auto& at_coord = curStep->getAtoms().coordinates;
    const auto& at_prop = curStep->getAtoms().properties;
    auto fmt = curStep->getFmt();
//    auto fmt_fun = curStep->getFormatter(fmt, AtomFmt::Bohr);
    float c, s, ic;
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
        cv = curStep->getCellVec() * curStep->getCellDim(AtomFmt::Angstrom);
        break;
    case AtomFmt::Bohr:
        break;
    }
    for(const Bond& bd:bonds){
        auto at_pos1 = at_coord[bd.at1];
        auto at_pos2 = at_coord[bd.at2];
        if (bd.diff[0]>0)     { at_pos2 += bd.diff[0]*cv[0]; }
        else if (bd.diff[0]<0){ at_pos1 -= bd.diff[0]*cv[0]; }
        if (bd.diff[1]>0)     { at_pos2 += bd.diff[1]*cv[1]; }
        else if (bd.diff[1]<0){ at_pos1 -= bd.diff[1]*cv[1]; }
        if (bd.diff[2]>0)     { at_pos2 += bd.diff[2]*cv[2]; }
        else if (bd.diff[2]<0){ at_pos1 -= bd.diff[2]*cv[2]; }
        auto bond_axis = at_pos1 - at_pos2;
        // FIXME: is this necessary?
//        if(fmt == AtomFmt::Crystal){
//            bond_axis = fmt_fun(bond_axis);
//        }
        auto bond_pos = (at_pos1+at_pos2)/2;
        const auto& col1 = bd.type ? bd.type->second : elements[bd.at1]->second.col;
        const auto& col2 = bd.type ? bd.type->second : elements[bd.at2]->second.col;
        // handle bonds parallel to x-axis
        if(std::abs(bond_axis[1])<std::numeric_limits<Vec::value_type>::epsilon()&&
           std::abs(bond_axis[2])<std::numeric_limits<Vec::value_type>::epsilon()){
            c = std::copysign(1.f, static_cast<float>(bond_axis[0]));
            bond_buffer.push_back({
                //mat3 with rotation and scaling
                {static_cast<float>(bd.dist)*c, 0., 0.,
                 0., bondRad, 0.,
                 0., 0., bondRad*c},
                //vec3 with position in modelspace
                {static_cast<float>(bond_pos[0]),
                 static_cast<float>(bond_pos[1]),
                 static_cast<float>(bond_pos[2])},
                //faux uvec4 with render criteria
                {static_cast<int16_t>(std::abs(bd.diff[0])),
                 static_cast<int16_t>(std::abs(bd.diff[1])),
                 static_cast<int16_t>(std::abs(bd.diff[2])),
                 static_cast<int16_t>(at_prop[bd.at1].flags[AtomFlag::Hidden] ||
                                      at_prop[bd.at2].flags[AtomFlag::Hidden])},
                 col1, col2});
        }else{
            // all other bonds
            auto rot_axis = -Vec_cross(bond_axis, x_axis);
            rot_axis /= Vec_length(rot_axis);
            c = static_cast<float>(Vec_dot(bond_axis, x_axis)/Vec_length(bond_axis));
            ic = 1-c;
            s = -std::sqrt(1-c*c);
            const auto dist = static_cast<float>(bd.dist);
            const float ax[3] = {static_cast<float>(rot_axis[0]),
                                 static_cast<float>(rot_axis[1]),
                                 static_cast<float>(rot_axis[2])};
            bond_buffer.push_back({
                //mat3 with rotation and scaling
                {dist*(ic*ax[0]*ax[0]+c),
                 dist*(ic*ax[0]*ax[1]-s*ax[2]),
                 dist*(ic*ax[0]*ax[2]+s*ax[1]),
                 bondRad*(ic*ax[1]*ax[0]+s*ax[2]),
                 bondRad*(ic*ax[1]*ax[1]+c),
                 bondRad*(ic*ax[1]*ax[2]-s*ax[0]),
                 bondRad*(ic*ax[2]*ax[0]-s*ax[1]),
                 bondRad*(ic*ax[2]*ax[1]+s*ax[0]),
                 bondRad*(ic*ax[2]*ax[2]+c)},
                //vec3 with position in modelspace
                {static_cast<float>(bond_pos[0]),
                 static_cast<float>(bond_pos[1]),
                 static_cast<float>(bond_pos[2])},
                //faux uvec4 with render criteria
                {static_cast<int16_t>(std::abs(bd.diff[0])),
                 static_cast<int16_t>(std::abs(bd.diff[1])),
                 static_cast<int16_t>(std::abs(bd.diff[2])),
                 static_cast<int16_t>(at_prop[bd.at1].flags[AtomFlag::Hidden] ||
                                      at_prop[bd.at2].flags[AtomFlag::Hidden])},
                //2*vec4 with colors
                col1, col2});
        }
    }
}
