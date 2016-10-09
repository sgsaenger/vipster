#include "glwidget.h"
#include <QKeyEvent>
#include "atom_model.h"
#include "bond_model.h"
#include <tuple>
#include <limits>
#include <cmath>
//DEBUGGING:
#include <iostream>
#include <chrono>

using namespace Vipster;

GLWidget::GLWidget(QWidget *parent):
    QOpenGLWidget(parent),
    cell_ibo(QOpenGLBuffer::IndexBuffer),
    aVo(true),
    bVo(true),
    cVo(true)
{
    QSurfaceFormat format;
    format.setVersion(3,3);
    format.setProfile(QSurfaceFormat::CoreProfile);
    this->setFormat(format);
}

GLWidget::~GLWidget()
{
    makeCurrent();

    vao.destroy();
    atom_vbo.destroy();
    for(auto vbo:bond_vbo) vbo.destroy();

    doneCurrent();
}

void GLWidget::initializeGL()
{
    initializeOpenGLFunctions();
    vao.create();
    if(vao.isCreated())vao.bind();
    atom_vbo.create();
    for(auto vbo:bond_vbo) vbo.create();

    cell_vbo.create();
    cell_ibo.create();
    ushort cell_idx[24] = {0,1,0,2,0,3,1,4,1,5,2,4,2,6,3,5,3,6,4,7,5,7,6,7};
    cell_ibo.bind();
    cell_ibo.allocate((void*)&cell_idx,sizeof(cell_idx));
    cell_ibo.release();

    sphere_vbo.create();
    sphere_vbo.bind();
    sphere_vbo.setUsagePattern(QOpenGLBuffer::StaticDraw);
    sphere_vbo.allocate((void*)&atom_model,4608*sizeof(float));
    sphere_vbo.release();

    torus_vbo.create();
    torus_vbo.bind();
    torus_vbo.setUsagePattern(QOpenGLBuffer::StaticDraw);
    torus_vbo.allocate((void*)&bond_model,288*sizeof(float));
    torus_vbo.release();

    glClearColor(1,1,1,1);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    glEnable(GL_BLEND);
    glEnable(GL_MULTISAMPLE);

    atom_shader.addShaderFromSourceFile(QOpenGLShader::Vertex,":/shaders/atom/atom.vert");
    atom_shader.addShaderFromSourceFile(QOpenGLShader::Fragment,":/shaders/atom/atom.frag");
    atom_shader.link();

    bond_shader.addShaderFromSourceFile(QOpenGLShader::Vertex,":/shaders/bond/bond.vert");
    bond_shader.addShaderFromSourceFile(QOpenGLShader::Fragment,":/shaders/bond/bond.frag");
    bond_shader.link();

    cell_shader.addShaderFromSourceFile(QOpenGLShader::Vertex,":/shaders/cell/cell.vert");
    cell_shader.addShaderFromSourceFile(QOpenGLShader::Fragment,":/shaders/cell/cell.frag");
    cell_shader.link();
}

void GLWidget::paintGL()
{
    vMatrix.setToIdentity();
    vMatrix.lookAt(QVector3D(0,0,10),QVector3D(),QVector3D(0,1,0));
    vMatrix.translate(0,0,0);
    drawAtoms();
    drawBonds();
    drawCell();
}

void GLWidget::drawAtoms(void)
{
    if(aVo){
        aVo=false;
        atom_vbo.bind();
        atom_vbo.allocate((void*)atom_buffer.data(),atom_buffer.size()*8*sizeof(float));
        atom_vbo.release();
    }
    atom_shader.bind();
    sphere_vbo.bind();
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,0,0);
    sphere_vbo.release();
    atom_shader.setUniformValue("vpMatrix",pMatrix*vMatrix*rMatrix);
    atom_shader.setUniformValue("rMatrix",rMatrix);
    atom_shader.setUniformValue("atom_fac",(GLfloat)0.5);

    atom_vbo.bind();
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,32,0);
    glVertexAttribDivisor(1,1);
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2,1,GL_FLOAT,GL_FALSE,32,(const GLvoid*)12);
    glVertexAttribDivisor(2,1);
    glEnableVertexAttribArray(3);
    glVertexAttribPointer(3,4,GL_FLOAT,GL_FALSE,32,(const GLvoid*)16);
    glVertexAttribDivisor(3,1);
    atom_vbo.release();
    glDrawArraysInstanced(GL_TRIANGLES,0,atom_model_npoly,atom_buffer.size());
    glDisableVertexAttribArray(0);
    glDisableVertexAttribArray(1);
    glVertexAttribDivisor(1,0);
    glDisableVertexAttribArray(2);
    glVertexAttribDivisor(2,0);
    glDisableVertexAttribArray(3);
    glVertexAttribDivisor(3,0);
    atom_shader.release();
}

void GLWidget::drawBonds(void)
{
    if(bVo){
        bVo=false;
        for(int i=0;i!=8;++i){
            bond_vbo[i].bind();
            bond_vbo[i].allocate((void*)bond_buffer[i].data(),bond_buffer[i].size()*24*sizeof(float));
            bond_vbo[i].release();
        }
    }
    bond_shader.bind();
    torus_vbo.bind();
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,0,0);
    torus_vbo.release();
    bond_shader.setUniformValue("vpMatrix",pMatrix*vMatrix*rMatrix);
    bond_shader.setUniformValue("rMatrix",rMatrix);

    bond_vbo[0].bind();
    for(long i=1;i!=7;++i){
       glEnableVertexAttribArray(i);
       glVertexAttribDivisor(i,1);
       glVertexAttribPointer(i,4,GL_FLOAT,GL_FALSE,24*sizeof(float),(const GLvoid*)((i-1)*16));
    }
    bond_vbo[0].release();
    glDrawArraysInstanced(GL_TRIANGLES,0,bond_model_npoly,bond_buffer[0].size());
    for(int i=1;i!=7;++i){
       glDisableVertexAttribArray(i);
       glVertexAttribDivisor(i,0);
    }
    bond_shader.release();
}

void GLWidget::drawCell(void)
{
    if(cVo){
        cVo=false;
        cell_vbo.bind();
        cell_vbo.allocate((void*)&cell_buffer,8*sizeof(Vec));
        cell_vbo.release();
    }
    cell_shader.bind();
    cell_shader.setUniformValue("vpMatrix",pMatrix*vMatrix*rMatrix);
    cell_vbo.bind();
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,0,(void*)0);
    cell_vbo.release();
    cell_ibo.bind();
    glDrawElements(GL_LINES,24,GL_UNSIGNED_SHORT,(void*)0);
    cell_ibo.release();
    glDisableVertexAttribArray(0);
    cell_shader.release();
}

void GLWidget::resizeGL(int w, int h)
{
    h==0?h=1:0;
    glViewport(0,0,w,h);
    float aspect = float(w)/h;
    pMatrix.setToIdentity();
    //pMatrix.perspective(60.0,aspect,0.001,1000);
    pMatrix.ortho(-10*aspect,10*aspect,-10,10,0,1000);
}

void GLWidget::setMode(int i,bool t)
{
    if(!t)return;
    i=i+1;
}

void GLWidget::setStep(const Step& step)
{
    aVo=bVo=cVo=true;
    //cell
    float cdm = step.getCellDim();
    Vec x,y,z,xy,xz,yz,xyz;
    x = step.getCellVec()[0]*cdm;
    y = step.getCellVec()[1]*cdm;
    z = step.getCellVec()[2]*cdm;
    xy = x+y;
    xz = x+z;
    yz = y+z;
    xyz = xy+z;
    cell_buffer = {{ Vec(),x,y,z,xy,xz,yz,xyz }};
    //atoms
    atom_buffer.reserve(step.getNat());
    atom_buffer.clear();
    for(const Atom& at:step.getAtoms()){
        PseEntry *pse = &step.pse[at.name];
        atom_buffer.push_back({{at.coord[0],at.coord[1],at.coord[2],pse->covr,
                          pse->col[0],pse->col[1],pse->col[2],pse->col[3]}});
    }
    //bonds
    Vec p1,p2,pos; //positions
    std::vector<float> c1,c2; //colors
    float c,s,ic; //cosine, sine, inverse cosine and angle for rot-mat
    float rad=1; //TODO: pull bond-radius from config
    Vec b_axis,r_axis; //bond, rotation axes
    Vec x_axis{{1,0,0}};
    for(int i=0;i!=8;++i){
        bond_buffer[i].reserve(step.getNat());
        bond_buffer[i].clear();
        for(const Bond& bd:step.getBonds()[i]){
            p1 = step.getAtom(bd.at1).coord;
            p2 = step.getAtom(bd.at2).coord;
            c1 = step.pse[step.getAtom(bd.at1).name].col;
            c2 = step.pse[step.getAtom(bd.at2).name].col;
            b_axis = p1-p2;
            pos = (p1+p2)/2;
            if(std::abs(b_axis[1])<std::numeric_limits<float>::epsilon()&&
                    std::abs(b_axis[2])<std::numeric_limits<float>::epsilon()){
                r_axis={{0,1,0}};
                c=std::copysign(1.,b_axis[0]);
                ic=1-c;
                s=0;
            }else{
                r_axis = -t_vec_cross(b_axis,x_axis);
                r_axis/=t_vec_length(r_axis);
                c = t_vec_dot(b_axis,x_axis)/t_vec_length(b_axis);
                ic=1-c;
                s = -std::sqrt(1-c*c);
            }
            bond_buffer[i].push_back({
                 bd.dist*(ic*r_axis[0]*r_axis[0]+c),
                 bd.dist*(ic*r_axis[0]*r_axis[1]-s*r_axis[2]),
                 bd.dist*(ic*r_axis[0]*r_axis[2]+s*r_axis[1]),
                 0,
                 rad*(ic*r_axis[1]*r_axis[0]+s*r_axis[2]),
                 rad*(ic*r_axis[1]*r_axis[1]+c),
                 rad*(ic*r_axis[1]*r_axis[2]-s*r_axis[0]),
                 0,
                 rad*(ic*r_axis[2]*r_axis[0]-s*r_axis[1]),
                 rad*(ic*r_axis[2]*r_axis[1]+s*r_axis[0]),
                 rad*(ic*r_axis[2]*r_axis[2]+c),
                 0,
                 pos[0],pos[1],pos[2],1,
                 c1[0],c1[1],c1[2],c1[3],
                 c2[0],c2[1],c2[2],c2[3]
             });
        }
    }
    update();
}

void GLWidget::keyPressEvent(QKeyEvent *e)
{
    switch(e->key()){
        case Qt::Key_Down:
            {QMatrix4x4 tmp;
            tmp.rotate(10,1,0,0);
            rMatrix = tmp*rMatrix;
            update();}
            break;
        case Qt::Key_Up:
            {QMatrix4x4 tmp;
            tmp.rotate(-10,1,0,0);
            rMatrix = tmp*rMatrix;
            update();}
            break;
        case Qt::Key_Left:
            {QMatrix4x4 tmp;
            tmp.rotate(-10,0,1,0);
            rMatrix = tmp*rMatrix;
            update();}
            break;
        case Qt::Key_Right:
            {QMatrix4x4 tmp;
            tmp.rotate(10,0,1,0);
            rMatrix = tmp*rMatrix;
            update();}
            break;
        default:
            break;
    }
}
