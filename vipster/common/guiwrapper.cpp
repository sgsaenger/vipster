#include <fstream>
#include <iostream>
#include <limits>

#include "guiwrapper.h"
#include "atom_model.h"
#include "bond_model.h"

using namespace Vipster;

#ifdef __EMSCRIPTEN__
std::string readShader(const std::string &filePath)
{
    std::string content;
    std::ifstream fileStream{filePath};

    if(!fileStream) throw std::invalid_argument{"Shader not found: "+filePath};

    content.assign(std::istreambuf_iterator<char>{fileStream},
                   std::istreambuf_iterator<char>{});
    return content;
}
#else
#include <QString>
#include <QFile>

std::string readShader(const std::string &filePath)
{
    QFile f(QString::fromStdString(filePath));
    f.open(QIODevice::ReadOnly);
    return f.readAll().toStdString();
}
#endif

void GuiWrapper::loadShader(GLuint &program, const std::string &header, std::string vertShaderStr, std::string fragShaderStr)
{
    GLuint vertShader = glCreateShader(GL_VERTEX_SHADER);
    GLuint fragShader = glCreateShader(GL_FRAGMENT_SHADER);
    GLint gl_ok = GL_FALSE;

    vertShaderStr = header + vertShaderStr;
    const char *vertShaderSrc = vertShaderStr.c_str();
    glShaderSource(vertShader, 1, &vertShaderSrc, nullptr);
    glCompileShader(vertShader);
    glGetShaderiv(vertShader, GL_COMPILE_STATUS, &gl_ok);
    if(gl_ok == 0){
        std::cout << "Vertex-Shader does not compile" << std::endl;
        GLint infoLen = 0;
        glGetShaderiv(vertShader, GL_INFO_LOG_LENGTH, &infoLen);
        std::vector<char> infoLog;
        infoLog.resize((infoLen > 1)?static_cast<size_t>(infoLen):1);
        glGetShaderInfoLog(vertShader, infoLen, nullptr, &infoLog[0]);
        std::cout << &infoLog[0] << std::endl;
        throw std::invalid_argument{"Vertex-Shader does not compile"};
    }

    fragShaderStr = header + fragShaderStr;
    const char *fragShaderSrc = fragShaderStr.c_str();
    glShaderSource(fragShader, 1, &fragShaderSrc, nullptr);
    glCompileShader(fragShader);
    glGetShaderiv(fragShader, GL_COMPILE_STATUS, &gl_ok);
    if(gl_ok == 0){
        std::cout << "Shader does not compile" << std::endl;
        GLint infoLen = 0;
        glGetShaderiv(fragShader, GL_INFO_LOG_LENGTH, &infoLen);
        std::vector<char> infoLog(infoLen > 1?static_cast<size_t>(infoLen):1);
        glGetShaderInfoLog(fragShader, infoLen, nullptr, &infoLog[0]);
        std::cout << &infoLog[0] << std::endl;
        throw std::invalid_argument{"Shader does not compile"};
    }

    program = glCreateProgram();
    glAttachShader(program, vertShader);
    glAttachShader(program, fragShader);
    glLinkProgram(program);
    glGetProgramiv(program, GL_LINK_STATUS, &gl_ok);
    if(gl_ok == 0){
        std::cout << "Program does not link" << std::endl;
        GLint infoLen = 0;
        std::vector<char> infoLog(infoLen > 1?static_cast<size_t>(infoLen):1);
        glGetProgramInfoLog(program, infoLen, nullptr, &infoLog[0]);
        std::cout << &infoLog[0] << std::endl;
        throw std::invalid_argument{"Program does not link"};
    }

    glDetachShader(program, vertShader);
    glDeleteShader(vertShader);
    glDetachShader(program, fragShader);
    glDeleteShader(fragShader);
}

void GuiWrapper::initGL(void)
{
    glClearColor(1,1,1,1);
#ifndef __EMSCRIPTEN__
    glEnable(GL_MULTISAMPLE);
#endif
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    //
    initAtomVAO();
    initBondVAO();
    initCellVAO();
    initSelVAO();
    initViewUBO();
    initViewMat();
}

void GuiWrapper::initShaders(const std::string& header, const std::string& folder)
{
    loadShader(atom_program, header,
               readShader(folder + "/atom.vert"),
               readShader(folder + "/atom.frag"));
    loadShader(bond_program, header,
               readShader(folder + "/bond.vert"),
               readShader(folder + "/bond.frag"));
    loadShader(cell_program, header,
               readShader(folder + "/cell.vert"),
               readShader(folder + "/cell.frag"));
    loadShader(sel_program, header,
               readShader(folder + "/select.vert"),
               readShader(folder + "/select.frag"));
}

void Vipster::guiMatScale(guiMat &m, float f)
{
    for(size_t i=0;i<4;i++){
        m[i*4+0]*=f;
        m[i*4+1]*=f;
        m[i*4+2]*=f;
    }
}

void Vipster::guiMatTranslate(guiMat &m, float x, float y, float z)
{
    //assuming 0 0 0 1 in last row of m
    m[3]+=x;
    m[7]+=y;
    m[11]+=z;
}

void Vipster::guiMatRot(guiMat &m, float a, float x, float y, float z)
{
    if(float_comp(a,0)){
        return;
    }
    float tmp = a * deg2rad;
    float s = std::sin(tmp);
    float c = std::cos(tmp);
    if(float_comp(x,0)){
        if(float_comp(y,0)){
            if(!float_comp(z,0)){
                // z-axis
                if (z<0){
                    s = -s;
                }
                m[0] = c * (tmp = m[0]) - s * m[4];
                m[4] = s * tmp + c * m[4];
                m[1] = c * (tmp = m[1]) - s * m[5];
                m[5] = s * tmp + c * m[5];
                m[2] = c * (tmp = m[2]) - s * m[6];
                m[6] = s * tmp + c * m[6];
                m[3] = c * (tmp = m[3]) - s * m[7];
                m[7] = s * tmp + c * m[7];
            }else{
                throw Error("guiMatRot: rotation axis-vector is zero");
            }
        }else if(float_comp(z,0)){
            // y-axis
            if (y<0) {
                s = -s;
            }
            m[0] = c * (tmp = m[0]) + s * m[8];
            m[8] = -s * tmp + c * m[8];
            m[1] = c * (tmp = m[1]) + s * m[9];
            m[9] = -s * tmp + c * m[9];
            m[2] = c * (tmp = m[2]) + s * m[10];
            m[10] = -s * tmp + c * m[10];
            m[3] = c * (tmp = m[3]) + s * m[11];
            m[11] = -s * tmp + c * m[11];
        }
    }else if(float_comp(y,0) && float_comp(z,0)){
        // x-axis
        if (x<0) {
            s = -s;
        }
        m[4] = c * (tmp = m[4]) - s * m[8];
        m[8] = s * tmp + c * m[8];
        m[5] = c * (tmp = m[5]) - s * m[9];
        m[9] = s * tmp + c * m[9];
        m[6] = c * (tmp = m[6]) - s * m[10];
        m[10] = s * tmp + c * m[10];
        m[7] = c * (tmp = m[7]) - s * m[11];
        m[11] = s * tmp + c * m[11];
    }else{
        // general rotation
        Vec axis{{x,y,z}};
        axis /= Vec_length(axis);
        Vec axismc = axis*(1-c);
        guiMat rotate{{c+axismc[0]*axis[0], axismc[1]*axis[0]-s*axis[2], axismc[2]*axis[0]+s*axis[1], 0,
                      axismc[0]*axis[1]+s*axis[2], c+axismc[1]*axis[1], axismc[2]*axis[1]-s*axis[0], 0,
                      axismc[0]*axis[2]-s*axis[1], axismc[1]*axis[2]+s*axis[0], c+axismc[2]*axis[2], 0,
                      0,0,0,1}};
        m = rotate * m;
    }
}

guiMat Vipster::guiMatMkOrtho(float l, float r, float b, float t, float n, float f)
{
    return guiMat{{2/(r-l), 0, 0, (r+l)/(l-r),
                   0, 2/(t-b), 0, (t+b)/(b-t),
                   0, 0, (2/(n-f)), ((f+n)/(n-f)),
                   0, 0, 0, 1}};
}

guiMat Vipster::guiMatMkLookAt(Vec eye, Vec target, Vec up)
{
    Vec dir = target - eye;
    dir /= Vec_length(dir);
    Vec r = Vec_cross(dir, up);
    r /= Vec_length(r);
    Vec u = Vec_cross(r, dir);
    return guiMat{{r[0], r[1], r[2], -Vec_dot(r, eye),
                  u[0], u[1], u[2], -Vec_dot(u, eye),
                  -dir[0], -dir[1], -dir[2], Vec_dot(dir, eye),
                  0, 0, 0, 1}};
}

guiMat Vipster::operator *=(guiMat &a, const guiMat &b)
{
    a = guiMat{{a[0]*b[0]+a[1]*b[4]+a[2]*b[ 8]+a[3]*b[12],
               a[0]*b[1]+a[1]*b[5]+a[2]*b[ 9]+a[3]*b[13],
               a[0]*b[2]+a[1]*b[6]+a[2]*b[10]+a[3]*b[14],
               a[0]*b[3]+a[1]*b[7]+a[2]*b[11]+a[3]*b[15],
               a[4]*b[0]+a[5]*b[4]+a[6]*b[ 8]+a[7]*b[12],
               a[4]*b[1]+a[5]*b[5]+a[6]*b[ 9]+a[7]*b[13],
               a[4]*b[2]+a[5]*b[6]+a[6]*b[10]+a[7]*b[14],
               a[4]*b[3]+a[5]*b[7]+a[6]*b[11]+a[7]*b[15],
               a[8]*b[0]+a[9]*b[4]+a[10]*b[ 8]+a[11]*b[12],
               a[8]*b[1]+a[9]*b[5]+a[10]*b[ 9]+a[11]*b[13],
               a[8]*b[2]+a[9]*b[6]+a[10]*b[10]+a[11]*b[14],
               a[8]*b[3]+a[9]*b[7]+a[10]*b[11]+a[11]*b[15],
               a[12]*b[0]+a[13]*b[4]+a[14]*b[ 8]+a[15]*b[12],
               a[12]*b[1]+a[13]*b[5]+a[14]*b[ 9]+a[15]*b[13],
               a[12]*b[2]+a[13]*b[6]+a[14]*b[10]+a[15]*b[14],
               a[12]*b[3]+a[13]*b[7]+a[14]*b[11]+a[15]*b[15]}};
    return a;
}

guiMat Vipster::operator *(guiMat a, const guiMat &b)
{
    return a*=b;
}

void GuiWrapper::initViewUBO(void)
{
    glGenBuffers(1, &view_ubo);
    glBindBuffer(GL_UNIFORM_BUFFER, view_ubo);
    glBufferData(GL_UNIFORM_BUFFER, 2*sizeof(guiMat), nullptr, GL_STATIC_DRAW);
    glBindBufferBase(GL_UNIFORM_BUFFER, 0, view_ubo);

    GLuint atom_loc = glGetUniformBlockIndex(atom_program, "viewMat");
    glUniformBlockBinding(atom_program, 0, atom_loc);
    GLuint bond_loc = glGetUniformBlockIndex(bond_program, "viewMat");
    glUniformBlockBinding(bond_program, 0, bond_loc);
    GLuint cell_loc = glGetUniformBlockIndex(cell_program, "viewMat");
    glUniformBlockBinding(cell_program, 0, cell_loc);
    GLuint sel_loc = glGetUniformBlockIndex(sel_program, "viewMat");
    glUniformBlockBinding(sel_program, 0, sel_loc);
}

void GuiWrapper::initAtomVAO(void)
{
    glGenVertexArrays(1, &atom_vao);
    glBindVertexArray(atom_vao);
    // upload and link sphere-mesh
    glGenBuffers(1, &sphere_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, sphere_vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(atom_model),
                 static_cast<const void*>(&atom_model), GL_STATIC_DRAW);
    auto vertexLoc = static_cast<GLuint>(glGetAttribLocation(atom_program, "vertex_modelspace"));
    glVertexAttribPointer(vertexLoc, 3, GL_FLOAT, GL_FALSE, 0, nullptr);
    glEnableVertexAttribArray(vertexLoc);
    // link atom positions
    glGenBuffers(1, &atom_pos_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, atom_pos_vbo);
    auto positionLoc = static_cast<GLuint>(glGetAttribLocation(atom_program, "position_modelspace"));
    glVertexAttribPointer(positionLoc, 3,
                          GL_FLOAT, GL_FALSE,
                          sizeof(Vec), nullptr);
    glVertexAttribDivisor(positionLoc, 1);
    glEnableVertexAttribArray(positionLoc);
    // link atom properties
    glGenBuffers(1, &atom_prop_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, atom_prop_vbo);
    auto scaleLoc = static_cast<GLuint>(glGetAttribLocation(atom_program, "scale_modelspace"));
    glVertexAttribPointer(scaleLoc, 1,
                          GL_FLOAT, GL_FALSE,
                          sizeof(atom_prop),
                          reinterpret_cast<const GLvoid*>(offsetof(atom_prop, rad)));
    glVertexAttribDivisor(scaleLoc, 1);
    glEnableVertexAttribArray(scaleLoc);
    auto colorLoc = static_cast<GLuint>(glGetAttribLocation(atom_program, "color_input"));
    glVertexAttribPointer(colorLoc, 4,
                          GL_UNSIGNED_BYTE, GL_TRUE,
                          sizeof(atom_prop),
                          reinterpret_cast<const GLvoid*>(offsetof(atom_prop, col)));
    glVertexAttribDivisor(colorLoc, 1);
    glEnableVertexAttribArray(colorLoc);
}

void GuiWrapper::initSelVAO(void)
{
    glGenVertexArrays(1, &sel_vao);
    glBindVertexArray(sel_vao);
    glBindBuffer(GL_ARRAY_BUFFER, sphere_vbo);
    auto vertexLoc = static_cast<GLuint>(glGetAttribLocation(atom_program, "vertex_modelspace"));
    glVertexAttribPointer(vertexLoc, 3, GL_FLOAT, GL_FALSE, 0, nullptr);
    glEnableVertexAttribArray(vertexLoc);
    glGenBuffers(1, &sel_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, sel_vbo);
    auto positionLoc = static_cast<GLuint>(glGetAttribLocation(atom_program, "position_modelspace"));
    glVertexAttribPointer(positionLoc, 3,
                          GL_FLOAT, GL_FALSE,
                          sizeof(sel_prop),
                          reinterpret_cast<const GLvoid*>(offsetof(sel_prop, pos)));
    glVertexAttribDivisor(positionLoc, 1);
    glEnableVertexAttribArray(positionLoc);
    auto scaleLoc = static_cast<GLuint>(glGetAttribLocation(atom_program, "scale_modelspace"));
    glVertexAttribPointer(scaleLoc, 1,
                          GL_FLOAT, GL_FALSE,
                          sizeof(sel_prop),
                          reinterpret_cast<const GLvoid*>(offsetof(sel_prop, rad)));
    glVertexAttribDivisor(scaleLoc, 1);
    glEnableVertexAttribArray(scaleLoc);
}

void GuiWrapper::initBondVAO(void)
{
    glGenVertexArrays(1, &bond_vao);
    glBindVertexArray(bond_vao);
    glGenBuffers(1, &torus_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, torus_vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(bond_model),
                 static_cast<const void*>(&bond_model), GL_STATIC_DRAW);
    auto vertexLoc = static_cast<GLuint>(glGetAttribLocation(bond_program, "vertex_modelspace"));
    glVertexAttribPointer(vertexLoc,3,GL_FLOAT,GL_FALSE,0,nullptr);
    glEnableVertexAttribArray(vertexLoc);
    glGenBuffers(1, &bond_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, bond_vbo);
    auto mMatLoc = static_cast<GLuint>(glGetAttribLocation(bond_program, "mMatrix"));
    glVertexAttribPointer(mMatLoc, 3,
                          GL_FLOAT,GL_FALSE,
                          sizeof(bond_prop), nullptr);
    glVertexAttribDivisor(mMatLoc, 1);
    glVertexAttribPointer(mMatLoc+1, 3,
                          GL_FLOAT,GL_FALSE,
                          sizeof(bond_prop),
                          reinterpret_cast<const GLvoid*>(offsetof(bond_prop, mat[3])));
    glVertexAttribDivisor(mMatLoc+1, 1);
    glVertexAttribPointer(mMatLoc+2, 3,
                          GL_FLOAT,GL_FALSE,
                          sizeof(bond_prop),
                          reinterpret_cast<const GLvoid*>(offsetof(bond_prop, mat[6])));
    glVertexAttribDivisor(mMatLoc+2, 1);
    glEnableVertexAttribArray(mMatLoc);
    glEnableVertexAttribArray(mMatLoc+1);
    glEnableVertexAttribArray(mMatLoc+2);
    auto positionLoc = static_cast<GLuint>(glGetAttribLocation(bond_program, "position_modelspace"));
    glVertexAttribPointer(positionLoc, 3,
                          GL_FLOAT,GL_FALSE,
                          sizeof(bond_prop),
                          reinterpret_cast<const GLvoid*>(offsetof(bond_prop, pos)));
    glVertexAttribDivisor(positionLoc, 1);
    glEnableVertexAttribArray(positionLoc);
    auto critLoc = static_cast<GLuint>(glGetAttribLocation(bond_program, "pbc_crit"));
    glVertexAttribIPointer(critLoc, 4,
                          GL_UNSIGNED_SHORT,
                          sizeof(bond_prop),
                          reinterpret_cast<const GLvoid*>(offsetof(bond_prop, mult)));
    glVertexAttribDivisor(critLoc, 1);
    glEnableVertexAttribArray(critLoc);
    auto color1Loc = static_cast<GLuint>(glGetAttribLocation(bond_program, "s1Color"));
    glVertexAttribPointer(color1Loc, 4,
                          GL_UNSIGNED_BYTE,GL_TRUE,
                          sizeof(bond_prop),
                          reinterpret_cast<const GLvoid*>(offsetof(bond_prop, col_a)));
    glVertexAttribDivisor(color1Loc, 1);
    glEnableVertexAttribArray(color1Loc);
    auto color2Loc = static_cast<GLuint>(glGetAttribLocation(bond_program, "s2Color"));
    glVertexAttribPointer(color2Loc, 4,
                          GL_UNSIGNED_BYTE,GL_TRUE,
                          sizeof(bond_prop),
                          reinterpret_cast<const GLvoid*>(offsetof(bond_prop, col_b)));
    glVertexAttribDivisor(color2Loc, 1);
    glEnableVertexAttribArray(color2Loc);
}

void GuiWrapper::initCellVAO(void)
{
    glGenVertexArrays(1, &cell_vao);
    glBindVertexArray(cell_vao);
    glGenBuffers(1, &cell_ibo);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cell_ibo);
    GLushort indices[24] = {0,1,0,2,0,3,1,4,1,5,2,4,2,6,3,5,3,6,4,7,5,7,6,7};
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), static_cast<void*>(indices), GL_STATIC_DRAW);
    glGenBuffers(1, &cell_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, cell_vbo);
    auto vertexLoc = static_cast<GLuint>(glGetAttribLocation(cell_program, "vertex_modelspace"));
    glVertexAttribPointer(vertexLoc, 3, GL_FLOAT, GL_FALSE, 0, nullptr);
    glEnableVertexAttribArray(vertexLoc);
}

void GuiWrapper::deleteGLObjects(void)
{
    glDeleteBuffers(1, &view_ubo);
    glDeleteProgram(atom_program);
    glDeleteBuffers(1, &atom_pos_vbo);
    glDeleteBuffers(1, &atom_prop_vbo);
    glDeleteBuffers(1, &sphere_vbo);
    glDeleteVertexArrays(1, &atom_vao);
    glDeleteProgram(bond_program);
    glDeleteBuffers(1, &bond_vbo);
    glDeleteBuffers(1, &torus_vbo);
    glDeleteVertexArrays(1, &bond_vao);
    glDeleteProgram(cell_program);
    glDeleteBuffers(1, &cell_vbo);
    glDeleteBuffers(1, &cell_ibo);
    glDeleteVertexArrays(1, &cell_vao);
}

void GuiWrapper::draw(void)
{
    glClearColor(1,1,1,1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    if(curStep->hasCell()){
        drawCell();
    }else{
        drawMol();
    }
}

void GuiWrapper::drawMol(void)
{
    Vec center = -curStep->getCenter(CdmFmt::Bohr);
    // atoms
    glBindVertexArray(atom_vao);
    glUseProgram(atom_program);
    auto offLocA = glGetUniformLocation(atom_program, "offset");
    auto facLoc = glGetUniformLocation(atom_program, "atom_fac");
    auto cellLocA = glGetUniformLocation(atom_program, "position_scale");
    auto toggleLoc = glGetUniformLocation(atom_program, "has_single_color");
    auto colorLoc = glGetUniformLocation(atom_program, "single_color");
    glUniform1f(facLoc, settings.atRadFac.val);
    glUniform3fv(offLocA, 1, center.data());
    glUniformMatrix3fv(cellLocA, 1, 0, cell_mat.data());
    glUniform1i(toggleLoc, 0);
    glDrawArraysInstanced(GL_TRIANGLES, 0,
                          atom_model_npoly,
                          static_cast<GLsizei>(atom_prop_buffer.size()));
    // bonds
    if(settings.showBonds.val){
        glBindVertexArray(bond_vao);
        glUseProgram(bond_program);
        auto offLocB = glGetUniformLocation(bond_program, "offset");
        auto pbcLoc = glGetUniformLocation(bond_program, "pbc_cell");
        auto multLoc = glGetUniformLocation(bond_program, "mult");
        auto cellLocB = glGetUniformLocation(bond_program, "position_scale");
        glUniform3ui(multLoc, 1, 1, 1);
        glUniformMatrix3fv(cellLocB, 1, 0, cell_mat.data());
        glUniform3fv(offLocB, 1, center.data());
        glUniform3ui(pbcLoc, 0, 0, 0);
        glDrawArraysInstanced(GL_TRIANGLES, 0,
                              bond_model_npoly,
                              static_cast<GLsizei>(bond_buffer.size()));
    }
    if(sel_buffer.size()){
        glBindVertexArray(sel_vao);
        glUseProgram(atom_program);
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

void GuiWrapper::drawCell(void)
{
    Vec off;
    Vec center = -curStep->getCenter(CdmFmt::Bohr);
    Mat cv = curStep->getCellVec() * curStep->getCellDim(CdmFmt::Bohr);
    center -= (mult[0]-1)*cv[0]/2.;
    center -= (mult[1]-1)*cv[1]/2.;
    center -= (mult[2]-1)*cv[2]/2.;
    // atoms
    glBindVertexArray(atom_vao);
    glUseProgram(atom_program);
    auto offLocA = glGetUniformLocation(atom_program, "offset");
    auto facLoc = glGetUniformLocation(atom_program, "atom_fac");
    auto cellLocA = glGetUniformLocation(atom_program, "position_scale");
    auto toggleLoc = glGetUniformLocation(atom_program, "has_single_color");
    auto colorLoc = glGetUniformLocation(atom_program, "single_color");
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
                                      static_cast<GLsizei>(atom_prop_buffer.size()));
            }
        }
    }
    // bonds
    if(settings.showBonds.val){
        glBindVertexArray(bond_vao);
        glUseProgram(bond_program);
        GLint offLocB = glGetUniformLocation(bond_program, "offset");
        GLint pbcLoc = glGetUniformLocation(bond_program, "pbc_cell");
        GLint multLoc = glGetUniformLocation(bond_program, "mult");
        GLint cellLocB = glGetUniformLocation(bond_program, "position_scale");
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
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cell_ibo);
        glUseProgram(cell_program);
        auto offLocC = glGetUniformLocation(cell_program, "offset");
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
    // selection
    if(sel_buffer.size()){
        glBindVertexArray(sel_vao);
        glUseProgram(atom_program);
        glUniform1f(facLoc, settings.atRadFac.val);
        glUniform3fv(offLocA, 1, (center).data());
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

void GuiWrapper::drawSel()
{
    // draw Atoms (as setup in atom_vao, shader locations must match)
    // with selection shader -> color by gl_InstanceID
    // to seperate Framebuffer
#ifndef __EMSCRIPTEN__
    glDisable(GL_MULTISAMPLE);
#endif
    glClearColor(1,1,1,0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    Vec center = -curStep->getCenter(CdmFmt::Bohr);
    glBindVertexArray(atom_vao);
    glUseProgram(sel_program);
    auto offLocA = glGetUniformLocation(sel_program, "offset");
    auto facLoc = glGetUniformLocation(sel_program, "atom_fac");
    auto cellLocA = glGetUniformLocation(sel_program, "position_scale");
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
                                          static_cast<GLsizei>(atom_prop_buffer.size()));
                }
            }
        }
    }else{
        glUniform3fv(offLocA, 1, center.data());
        glDrawArraysInstanced(GL_TRIANGLES, 0,
                              atom_model_npoly,
                              static_cast<GLsizei>(atom_prop_buffer.size()));
    }
#ifndef __EMSCRIPTEN__
    glEnable(GL_MULTISAMPLE);
#endif
}

void GuiWrapper::updateStepBuffers(StepProper* step, bool draw_bonds)
{
    if (step!=nullptr) {
        curStep = step;
    }
    //cell
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

    //atoms
    atom_prop_buffer.clear();
    atom_prop_buffer.reserve(curStep->getNat());
    for (const PseEntry* at:curStep->getPseEntries()){
        atom_prop_buffer.push_back({at->covr, at->col});
    }
    atoms_changed = true;

    //bonds
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

void GuiWrapper::updateSelBuffers(StepSelection* sel)
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
}

void GuiWrapper::updateVBOs(void)
{
    if (atoms_changed) {
        glBindBuffer(GL_ARRAY_BUFFER, atom_pos_vbo);
        auto nat = atom_prop_buffer.size();
        if (nat != 0u) {
            glBufferData(GL_ARRAY_BUFFER, static_cast<GLsizeiptr>(nat*sizeof(Vec)),
                         static_cast<const void*>(curStep->getCoords().data()), GL_STREAM_DRAW);
        } else {
            glBufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_STREAM_DRAW);
        }
        glBindBuffer(GL_ARRAY_BUFFER, atom_prop_vbo);
        if (nat != 0u) {
            glBufferData(GL_ARRAY_BUFFER, static_cast<GLsizeiptr>(nat*sizeof(atom_prop)),
                         static_cast<void*>(atom_prop_buffer.data()), GL_STREAM_DRAW);
        } else {
            glBufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_STREAM_DRAW);
        }
        atoms_changed = false;
    }
    if (sel_changed) {
        glBindBuffer(GL_ARRAY_BUFFER, sel_vbo);
        if(!sel_buffer.empty()){
            glBufferData(GL_ARRAY_BUFFER, static_cast<GLsizeiptr>(sel_buffer.size()*sizeof(sel_prop)),
                         static_cast<const void*>(sel_buffer.data()), GL_STREAM_DRAW);
        }else{
            glBufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_STREAM_DRAW);
        }
        sel_changed = false;
    }
    if (bonds_changed) {
        glBindBuffer(GL_ARRAY_BUFFER, bond_vbo);
        if (!bond_buffer.empty()) {
            glBufferData(GL_ARRAY_BUFFER, static_cast<GLsizeiptr>(bond_buffer.size()*sizeof(bond_prop)),
                         static_cast<void*>(bond_buffer.data()), GL_STREAM_DRAW);
        } else {
            glBufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_STREAM_DRAW);
        }
        bonds_changed = false;
    }
    if (cell_changed) {
        glBindBuffer(GL_ARRAY_BUFFER, cell_vbo);
        glBufferData(GL_ARRAY_BUFFER, 8*sizeof(Vec),
                     static_cast<void*>(cell_buffer.data()), GL_STREAM_DRAW);
        cell_changed = false;
    }
}

void GuiWrapper::updateViewUBO(void)
{
    if(rMatChanged){
        glBindBuffer(GL_UNIFORM_BUFFER, view_ubo);
        glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(guiMat),
                        (pMat*vMat*rMat).data());
        glBufferSubData(GL_UNIFORM_BUFFER, sizeof(guiMat), sizeof(guiMat), rMat.data());
        rMatChanged = vMatChanged = pMatChanged = false;
    }else if (pMatChanged || vMatChanged){
        glBindBuffer(GL_UNIFORM_BUFFER, view_ubo);
        glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(guiMat),
                        (pMat*vMat*rMat).data());
        vMatChanged = pMatChanged = false;
    }
}

void GuiWrapper::initViewMat(void)
{
    pMat = guiMatMkOrtho(-15, 15, -10, 10, -100, 1000);
    rMat = {{1,0,0,0,
             0,1,0,0,
             0,0,1,0,
             0,0,0,1}};
    vMat = guiMatMkLookAt({{0,0,10}}, {{0,0,0}}, {{0,1,0}});
    pMatChanged = rMatChanged = vMatChanged = true;
}

void GuiWrapper::resizeViewMat(int w, int h)
{
    h==0?h=1:0;
    glViewport(0,0,w,h);
    float aspect = float(w)/h;
    pMat = guiMatMkOrtho(-10*aspect, 10*aspect, -10, 10, -100, 1000);
    pMatChanged = true;
}

void GuiWrapper::zoomViewMat(int i)
{
    guiMatScale(vMat, i>0?1.1f:0.9f);
    vMatChanged = true;
}

void GuiWrapper::rotateViewMat(float x, float y, float z)
{
    guiMatRot(rMat, x, 0, 1, 0);
    guiMatRot(rMat, y, 1, 0, 0);
    guiMatRot(rMat, z, 0, 0, 1);
    rMatChanged = true;
}

void GuiWrapper::translateViewMat(float x, float y, float z)
{
    guiMatTranslate(vMat, x/10.f, y/10.f, z/10.f);
    vMatChanged = true;
}

void GuiWrapper::alignViewMat(alignDir d)
{
    switch (d) {
    case alignDir::x:
        rMat = {{ 0, 1, 0, 0,
                  0, 0, 1, 0,
                  1, 0, 0, 0,
                  0, 0, 0, 1}};
        break;
    case alignDir::mx:
        rMat = {{ 0,-1, 0, 0,
                  0, 0, 1, 0,
                 -1, 0, 0, 0,
                  0, 0, 0, 1}};
        break;
    case alignDir::y:
        rMat = {{-1, 0, 0, 0,
                  0, 0, 1, 0,
                  0, 1, 0, 0,
                  0, 0, 0, 1}};
        break;
    case alignDir::my:
        rMat = {{ 1, 0, 0, 0,
                  0, 0, 1, 0,
                  0,-1, 0, 0,
                  0, 0, 0, 1}};
        break;
    case alignDir::z:
        rMat = {{ 1, 0, 0, 0,
                  0, 1, 0, 0,
                  0, 0, 1, 0,
                  0, 0, 0, 1}};
        break;
    case alignDir::mz:
        rMat = {{-1, 0, 0, 0,
                  0, 1, 0, 0,
                  0, 0,-1, 0,
                  0, 0, 0, 1}};
        break;
    }
    rMatChanged = true;
}

Mat GuiWrapper::getAxes()
{
    Mat tmp;
    tmp[0] =  Vec{rMat[0], rMat[1], rMat[2]};
    tmp[1] = -Vec{rMat[4], rMat[5], rMat[6]};
    tmp[2] =  Vec{rMat[8], rMat[9], rMat[10]};
    return tmp;
}
