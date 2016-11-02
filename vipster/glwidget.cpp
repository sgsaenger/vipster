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
    QOpenGLWidget(parent)
{
    QSurfaceFormat format;
    format.setVersion(3,3);
    format.setSamples(8);
    format.setAlphaBufferSize(8);
    format.setProfile(QSurfaceFormat::CoreProfile);
    this->setFormat(format);
}

GLWidget::~GLWidget()
{
    makeCurrent();

    vao.destroy();
    atom_vbo.destroy();
    bond_vbo.destroy();

    doneCurrent();
}

void GLWidget::initializeGL()
{
    initializeOpenGLFunctions();
    glEnable(GL_MULTISAMPLE);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    vao.create();
    if(vao.isCreated())vao.bind();
    atom_vbo.create();
    bond_vbo.create();

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
    vMatrix.translate(xshift, yshift, 0);
    vMatrix.scale(distance);
    drawAtoms();
    drawBonds();
    drawCell();
}

void GLWidget::drawAtoms()
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
    glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,8*sizeof(float),0);
    glVertexAttribDivisor(1,1);
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2,1,GL_FLOAT,GL_FALSE,8*sizeof(float),(const GLvoid*)(3*sizeof(float)));
    glVertexAttribDivisor(2,1);
    glEnableVertexAttribArray(3);
    glVertexAttribPointer(3,4,GL_FLOAT,GL_FALSE,8*sizeof(float),(const GLvoid*)(4*sizeof(float)));
    glVertexAttribDivisor(3,1);
    atom_vbo.release();
    Vec center = curStep->getCenter();
    Vec xvec = curStep->getCellVec()[0]*curStep->getCellDim();
    Vec yvec = curStep->getCellVec()[1]*curStep->getCellDim();
    Vec zvec = curStep->getCellVec()[2]*curStep->getCellDim();
    center += (mult[0]-1)*xvec/2.;
    center += (mult[1]-1)*yvec/2.;
    center += (mult[2]-1)*zvec/2.;
    for(int x=0;x!=mult[0];++x){
        for(int y=0;y!=mult[1];++y){
            for(int z=0;z!=mult[2];++z){
                Vec off = -center + x*xvec + y*yvec + z*zvec;
                atom_shader.setUniformValue("offset", QVector3D(off[0], off[1], off[2]));
                glDrawArraysInstanced(GL_TRIANGLES,0,atom_model_npoly,atom_buffer.size());
            }
        }
    }
    glDisableVertexAttribArray(0);
    glDisableVertexAttribArray(1);
    glVertexAttribDivisor(1,0);
    glDisableVertexAttribArray(2);
    glVertexAttribDivisor(2,0);
    glDisableVertexAttribArray(3);
    glVertexAttribDivisor(3,0);
    atom_shader.release();
}

void GLWidget::drawBonds()
{
    if(bVo){
        bVo=false;
        bond_vbo.bind();
        bond_vbo.allocate((void*)bond_buffer.data(),bond_buffer.size()*24*sizeof(float));
        bond_vbo.release();
    }
    bond_shader.bind();
    torus_vbo.bind();
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,0,0);
    torus_vbo.release();
    bond_shader.setUniformValue("vpMatrix",pMatrix*vMatrix*rMatrix);
    bond_shader.setUniformValue("rMatrix",rMatrix);
    bond_shader.setUniformValue("mult",mult[0], mult[1], mult[2]);

    bond_vbo.bind();
    for(long i=1;i!=8;++i){
        glEnableVertexAttribArray(i);
        glVertexAttribDivisor(i,1);
    }
    glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,24*sizeof(float),0);
    glVertexAttribPointer(2,3,GL_FLOAT,GL_FALSE,24*sizeof(float),(const GLvoid*)(3*sizeof(float)));
    glVertexAttribPointer(3,3,GL_FLOAT,GL_FALSE,24*sizeof(float),(const GLvoid*)(6*sizeof(float)));
    glVertexAttribPointer(4,3,GL_FLOAT,GL_FALSE,24*sizeof(float),(const GLvoid*)(9*sizeof(float)));
    glVertexAttribPointer(5,4,GL_FLOAT,GL_FALSE,24*sizeof(float),(const GLvoid*)(12*sizeof(float)));
    glVertexAttribPointer(6,4,GL_FLOAT,GL_FALSE,24*sizeof(float),(const GLvoid*)(16*sizeof(float)));
    glVertexAttribPointer(7,4,GL_FLOAT,GL_FALSE,24*sizeof(float),(const GLvoid*)(20*sizeof(float)));
    bond_vbo.release();
    Vec center = curStep->getCenter();
    Vec xvec = curStep->getCellVec()[0]*curStep->getCellDim();
    Vec yvec = curStep->getCellVec()[1]*curStep->getCellDim();
    Vec zvec = curStep->getCellVec()[2]*curStep->getCellDim();
    center += (mult[0]-1)*xvec/2.;
    center += (mult[1]-1)*yvec/2.;
    center += (mult[2]-1)*zvec/2.;
    for(int x=0;x!=mult[0];++x){
        for(int y=0;y!=mult[1];++y){
            for(int z=0;z!=mult[2];++z){
                Vec off = -center + x*xvec + y*yvec + z*zvec;
                bond_shader.setUniformValue("offset", QVector3D(off[0], off[1], off[2]));
                bond_shader.setUniformValue("pbc_cell", QVector3D(x,y,z));
                glDrawArraysInstanced(GL_TRIANGLES,0,bond_model_npoly,bond_buffer.size());
            }
        }
    }
    for(int i=1;i!=8;++i){
       glDisableVertexAttribArray(i);
       glVertexAttribDivisor(i,0);
    }
    bond_shader.release();
}

void GLWidget::drawCell()
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
    Vec center = curStep->getCenter();
    Vec xvec = curStep->getCellVec()[0]*curStep->getCellDim();
    Vec yvec = curStep->getCellVec()[1]*curStep->getCellDim();
    Vec zvec = curStep->getCellVec()[2]*curStep->getCellDim();
    center += (mult[0]-1)*xvec/2.;
    center += (mult[1]-1)*yvec/2.;
    center += (mult[2]-1)*zvec/2.;
    for(int x=0;x!=mult[0];++x){
        for(int y=0;y!=mult[1];++y){
            for(int z=0;z!=mult[2];++z){
                Vec off = -center + x*xvec + y*yvec + z*zvec;
                cell_shader.setUniformValue("offset", QVector3D(off[0], off[1], off[2]));
                glDrawElements(GL_LINES,24,GL_UNSIGNED_SHORT,(void*)0);
            }
        }
    }
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
}

void GLWidget::setMult(int i)
{
    if(QObject::sender()->objectName() == "xMultBox"){ mult[0] = i; }
    else if(QObject::sender()->objectName() == "yMultBox"){ mult[1] = i; }
    else if(QObject::sender()->objectName() == "zMultBox"){ mult[2] = i; }
    update();
}

void GLWidget::setStep(const Step& step)
{
    aVo=bVo=cVo=true;
    curStep = &step;
    //cell
    float cdm = step.getCellDim();
    Vec x_vec,y_vec,z_vec;
    x_vec = step.getCellVec()[0]*cdm;
    y_vec = step.getCellVec()[1]*cdm;
    z_vec = step.getCellVec()[2]*cdm;
    cell_buffer = {{ Vec(),x_vec,y_vec,z_vec,x_vec+y_vec,x_vec+z_vec,
                     y_vec+z_vec,x_vec+y_vec+z_vec }};
    //atoms
    const std::vector<Atom>& atoms = step.getAtoms();
    atom_buffer.reserve(step.getNat());
    atom_buffer.clear();
    for(const Atom& at:atoms){
        PseEntry *pse = &step.pse[at.name];
        atom_buffer.push_back({{at.coord[0],at.coord[1],at.coord[2],pse->covr,
                          pse->col[0],pse->col[1],pse->col[2],pse->col[3]}});
    }
    //bonds
    Vec p1, p2, pos; //positions
    std::vector<float> c1, c2; //colors
    float c, s, ic; //cosine, sine, inverse cosine and angle for rot-mat
    float rad = 0.53; //TODO: pull bond-radius from config
    Vec b_axis, r_axis; //bond, rotation axes
    constexpr Vec x_axis{{1,0,0}};
    bond_buffer.reserve(step.getNat());
    bond_buffer.clear();
    for(const Bond& bd:step.getBondsCell()){
        const Atom &at1 = atoms[bd.at1];
        const Atom &at2 = atoms[bd.at2];
        c1 = step.pse[at1.name].col;
        std::vector<float> &c1 = step.pse[at1.name].col;
        std::vector<float> &c2 = step.pse[at2.name].col;
        p1 = at1.coord;
        p2 = at2.coord;
        if (bd.xdiff>0){ p2 += bd.xdiff * x_vec; }else if(bd.xdiff<0){ p1 -= bd.xdiff * x_vec; }
        if (bd.ydiff>0){ p2 += bd.ydiff * y_vec; }else if(bd.ydiff<0){ p1 -= bd.ydiff * y_vec; }
        if (bd.zdiff>0){ p2 += bd.zdiff * z_vec; }else if(bd.zdiff<0){ p1 -= bd.zdiff * z_vec; }
        b_axis = p1-p2;
        pos = (p1+p2)/2;
        if(std::abs(b_axis[1])<std::numeric_limits<float>::epsilon()&&
                std::abs(b_axis[2])<std::numeric_limits<float>::epsilon()){
            r_axis = {{0,1,0}};
            c = std::copysign(1.,b_axis[0]);
            ic = 1-c;
            s = 0;
        }else{
            r_axis = -Vec_cross(b_axis,x_axis);
            r_axis /= Vec_length(r_axis);
            c = Vec_dot(b_axis,x_axis)/Vec_length(b_axis);
            ic = 1-c;
            s = -std::sqrt(1-c*c);
        }
        bond_buffer.push_back({
             //mat3 with rotation and scaling
             bd.dist*(ic*r_axis[0]*r_axis[0]+c),
             bd.dist*(ic*r_axis[0]*r_axis[1]-s*r_axis[2]),
             bd.dist*(ic*r_axis[0]*r_axis[2]+s*r_axis[1]),
             rad*(ic*r_axis[1]*r_axis[0]+s*r_axis[2]),
             rad*(ic*r_axis[1]*r_axis[1]+c),
             rad*(ic*r_axis[1]*r_axis[2]-s*r_axis[0]),
             rad*(ic*r_axis[2]*r_axis[0]-s*r_axis[1]),
             rad*(ic*r_axis[2]*r_axis[1]+s*r_axis[0]),
             rad*(ic*r_axis[2]*r_axis[2]+c),
             //vec3 with position in modelspace
             pos[0],pos[1],pos[2],
             //faux vec4 with integral pbc information
             (float)std::abs(bd.xdiff),
             (float)std::abs(bd.ydiff),
             (float)std::abs(bd.zdiff),
             //padding float that tells if non-pbc bond
             (float)!(bd.xdiff||bd.ydiff||bd.zdiff),
             //2*vec4 with colors
             c1[0],c1[1],c1[2],c1[3],
             c2[0],c2[1],c2[2],c2[3]
        });
    }
    update();
}

void GLWidget::setCamera(int i)
{
    switch(i){
    case -2: // +x
        rMatrix = QMatrix4x4{ 0, 1, 0, 0,
                              0, 0, 1, 0,
                              1, 0, 0, 0,
                              0, 0, 0, 1};
        break;
    case -5: // -x
        rMatrix = QMatrix4x4{ 0,-1, 0, 0,
                              0, 0, 1, 0,
                             -1, 0, 0, 0,
                              0, 0, 0, 1};
        break;
    case -3: // +y
        rMatrix = QMatrix4x4{-1, 0, 0, 0,
                              0, 0, 1, 0,
                              0, 1, 0, 0,
                              0, 0, 0, 1};
        break;
    case -6: // -y
        rMatrix = QMatrix4x4{ 1, 0, 0, 0,
                              0, 0, 1, 0,
                              0,-1, 0, 0,
                              0, 0, 0, 1};
        break;
    case -4: // +z
        rMatrix.setToIdentity();
        break;
    case -7: // -z
        rMatrix = QMatrix4x4{-1, 0, 0, 0,
                              0, 1, 0, 0,
                              0, 0,-1, 0,
                              0, 0, 0, 1};
        break;
    }
    std::cout << i << std::endl;
    update();
}

void GLWidget::keyPressEvent(QKeyEvent *e)
{
    switch(e->key()){
        case Qt::Key_Down:
            {QMatrix4x4 tmp;
            tmp.rotate(10,1,0,0);
            rMatrix = tmp*rMatrix;
            update();
            e->accept();}
            break;
        case Qt::Key_Up:
            {QMatrix4x4 tmp;
            tmp.rotate(-10,1,0,0);
            rMatrix = tmp*rMatrix;
            update();
            e->accept();}
            break;
        case Qt::Key_Left:
            {QMatrix4x4 tmp;
            tmp.rotate(-10,0,1,0);
            rMatrix = tmp*rMatrix;
            update();
            e->accept();}
            break;
        case Qt::Key_Right:
            {QMatrix4x4 tmp;
            tmp.rotate(10,0,1,0);
            rMatrix = tmp*rMatrix;
            update();
            e->accept();}
            break;
        default:
            break;
    }
}

void GLWidget::wheelEvent(QWheelEvent *e)
{
    int delta = e->angleDelta().y();
    distance *= delta>0 ? 1.1 : 0.9;
    update();
    e->accept();
}

void GLWidget::mousePressEvent(QMouseEvent *e)
{
    setFocus();
    if (!(e->buttons()&7)) return;
    mousePos = e->pos();
    switch(mouseMode){
        case MouseMode::Camera:
            if(e->button() == Qt::MouseButton::RightButton){
                rMatrix.setToIdentity();
                update();
            }
            break;
        case MouseMode::Select:
            break;
        case MouseMode::Modify:
            break;
    }
}

void GLWidget::mouseMoveEvent(QMouseEvent *e)
{
    if (!(e->buttons()&7)) return;
    QPoint delta = e->pos() - mousePos;
    switch(mouseMode){
        case MouseMode::Camera:
            if(e->buttons() & Qt::MouseButton::LeftButton){
                QMatrix4x4 tmp;
                tmp.rotate(delta.x(),0,1,0);
                tmp.rotate(delta.y(),1,0,0);
                rMatrix = tmp * rMatrix;
                update();
            }else if(e->buttons() & Qt::MouseButton::MiddleButton){
                xshift += delta.x() / 10.;
                yshift -= delta.y() / 10.;
                update();
            }
            break;
        case MouseMode::Select:
            break;
        case MouseMode::Modify:
            break;
    }
    mousePos = e->pos();
    e->accept();
}

void GLWidget::mouseReleaseEvent(QMouseEvent *e)
{
    if (!(e->buttons()&7)) return;
    switch(mouseMode){
        case MouseMode::Camera:
            break;
        case MouseMode::Select:
            break;
        case MouseMode::Modify:
            break;
    }
}
