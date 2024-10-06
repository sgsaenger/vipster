#include "seldata.h"
#include "atom_model.h"

using namespace Vipster;

decltype(GUI::SelData::shader_map) GUI::SelData::shader_map;

GUI::SelData::SelData(Step::const_selection const* sel)
    : curSel{sel}
{}

GUI::SelData::SelData(SelData&& dat)
    : Data{std::move(dat)}
//      sel_buffer{std::move(dat.sel_buffer)},
//      cell_mat{dat.cell_mat},
//      curSel{dat.curSel}
{
    std::swap(sel_buffer, dat.sel_buffer);
    std::swap(cell_mat, dat.cell_mat);
    std::swap(curSel, dat.curSel);
    std::swap(atRadFac, dat.atRadFac);
    std::swap(color, dat.color);
    std::swap(object_map, object_map);
}

void GUI::SelData::initGL(void *context)
{
    auto &shader = shader_map[context];
    if(!shader.initialized){
        initShader(global_map[context], shader);
        shader.initialized = true;
    }
    auto &objects = object_map[context];
    if(!objects.initialized){
        glGenVertexArrays(1, &objects.vao);
        glGenBuffers(1, &objects.vbo);
        initVAO(global_map[context], objects, shader);
        objects.initialized = true;
    }
    glBindVertexArray(0);
}

void GUI::SelData::initShader(GlobalContext& globals, shader &shader)
{
    shader.program = loadShader(globals, "/selection.vert", "/selection.frag");
    READATTRIB(shader, vertex)
    READATTRIB(shader, position)
    READATTRIB(shader, vert_scale)
    READATTRIB(shader, pbc_crit)
    READUNIFORM(shader, color)
    READUNIFORM(shader, offset)
    READUNIFORM(shader, pos_scale)
    READUNIFORM(shader, scale_fac)
    READUNIFORM(shader, mult)
}

void GUI::SelData::initVAO(GlobalContext& globals, ObjectContext &objects, shader &shader)
{
    glBindVertexArray(objects.vao);

    glBindBuffer(GL_ARRAY_BUFFER, globals.sphere_vbo);
    glVertexAttribPointer(shader.vertex, 3, GL_FLOAT, GL_FALSE, 0, nullptr);
    glEnableVertexAttribArray(shader.vertex);

    glBindBuffer(GL_ARRAY_BUFFER, objects.vbo);

    // ATOM POSITIONS
    glVertexAttribPointer(shader.position, 3,
                          GL_FLOAT, GL_FALSE,
                          sizeof(SelProp),
                          reinterpret_cast<const GLvoid*>(offsetof(SelProp, pos)));
    glVertexAttribDivisor(shader.position, 1);
    glEnableVertexAttribArray(shader.position);

    // ATOM PROPERTIES
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
}

GUI::SelData::~SelData()
{
    for(auto &[context, objects]: object_map){
        if(!objects.initialized) continue;
        glDeleteVertexArrays(1, &objects.vao);
        glDeleteBuffers(1, &objects.vbo);
    }
}

void GUI::SelData::draw(const Vec &off, const PBCVec &mult,
                        const Mat &, bool, void *context)
{
    if(sel_buffer.size()){
        glBindVertexArray(object_map[context].vao);
        auto &shader = shader_map[context];
        glUseProgram(shader.program);
        glUniform1f(shader.scale_fac, atRadFac);
        glUniform3f(shader.offset,
                    static_cast<float>(off[0]),
                    static_cast<float>(off[1]),
                    static_cast<float>(off[2]));
        glUniformMatrix3fv(shader.pos_scale, 1, 0, cell_mat.data());
        glUniform4ui(shader.color, color[0], color[1], color[2], color[3]);
        glUniform3i(shader.mult, mult[0], mult[1], mult[2]);
        glDrawArraysInstanced(GL_TRIANGLES, 0,
                              atom_model_npoly,
                              static_cast<GLsizei>(sel_buffer.size()));
    }
}

void GUI::SelData::updateGL(void *context)
{
    glBindBuffer(GL_ARRAY_BUFFER, object_map[context].vbo);
    if(!sel_buffer.empty()){
        glBufferData(GL_ARRAY_BUFFER, static_cast<GLsizeiptr>(sel_buffer.size()*sizeof(SelProp)),
                     static_cast<const void*>(sel_buffer.data()), GL_STREAM_DRAW);
    }else{
        glBufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_STREAM_DRAW);
    }
}

void GUI::SelData::update(Step::const_selection const* sel, bool useVdW, float atRadFac)
{
    this->atRadFac = atRadFac;
    curSel = sel;
    sel_buffer.clear();
    if(!curSel){
        return;
    }
    sel_buffer.reserve(curSel->getNat()); // too small, but better than nothing
    if(useVdW){
        auto it = curSel->cbegin();
        while(it != curSel->cend()){
            const Vec& pos = it->coord;
            const SizeVec& off = it->off;
            sel_buffer.push_back({{static_cast<float>(pos[0]),
                                   static_cast<float>(pos[1]),
                                   static_cast<float>(pos[2])},
                                  static_cast<float>(it->type->vdwr*1.3),
                                  {static_cast<int16_t>(off[0]),
                                   static_cast<int16_t>(off[1]),
                                   static_cast<int16_t>(off[2])},
                                 });
            ++it;
        }
    }else{
        auto it = curSel->cbegin();
        while(it != curSel->cend()){
            const Vec& pos = it->coord;
            const SizeVec& off = it->off;
            sel_buffer.push_back({{static_cast<float>(pos[0]),
                                   static_cast<float>(pos[1]),
                                   static_cast<float>(pos[2])},
                                  static_cast<float>(it->type->covr*1.3),
                                  {static_cast<int16_t>(off[0]),
                                   static_cast<int16_t>(off[1]),
                                   static_cast<int16_t>(off[2])},
                                 });
            ++it;
        }
    }
    for(auto &[context, state]: instance_map){
        state.synchronized = false;
    }

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
        tmp_mat *= curSel->getCellDim(AtomFmt::Bohr);
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
}
