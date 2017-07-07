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

    doneCurrent();
}

void GLWidget::initializeGL()
{
    std::cout << "hallo!" << std::endl;
    initializeOpenGLFunctions();
    glClearColor(1,1,1,1);
    glEnable(GL_MULTISAMPLE);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    loadShader(atom_program, "# version 330\n",
                   readShader(":/shaders/atom/atom.vert"),
                   readShader(":/shaders/atom/atom.frag"));
    glUseProgram(atom_program);
    //TODO: remove temps
    glUniform1f(glGetUniformLocation(atom_program, "atom_fac"), 0.5);
    //
    initAtomVAO();
    //
    rMat = {{1,0,0,0,
                0,1,0,0,
                0,0,1,0,
                0,0,0,1}};
    vMat = guiMatMkLookAt({{0,0,10}}, {{0,0,0}}, {{0,1,0}});
    pMatChanged = rMatChanged = vMatChanged = true;
}

void GLWidget::paintGL()
{
    glUseProgram(atom_program);
    if(rMatChanged){
        glUniformMatrix4fv(glGetUniformLocation(atom_program, "rMatrix"),
                           1, true, rMat.data());
        glUniformMatrix4fv(glGetUniformLocation(atom_program, "vpMatrix"),
                           1, true, (pMat*vMat*rMat).data());
        rMatChanged = vMatChanged = pMatChanged = false;
    }else if (pMatChanged || vMatChanged){
        glUniformMatrix4fv(glGetUniformLocation(atom_program, "vpMatrix"),
                           1, true, (pMat*vMat*rMat).data());
        vMatChanged = pMatChanged = false;
    }
    if(atoms_changed){
        glBindBuffer(GL_ARRAY_BUFFER, atom_vbo);
        glBufferData(GL_ARRAY_BUFFER, atom_buffer.size()*8*sizeof(float),
                     (void*)atom_buffer.data(), GL_STREAM_DRAW);
        atoms_changed=false;
    }
    //
    Vec center = curStep->getCenter();
    Vec xvec = curStep->getCellVec()[0]*curStep->getCellDim();
    Vec yvec = curStep->getCellVec()[1]*curStep->getCellDim();
    Vec zvec = curStep->getCellVec()[2]*curStep->getCellDim();
    center += (mult[0]-1)*xvec/2.;
    center += (mult[1]-1)*yvec/2.;
    center += (mult[2]-1)*zvec/2.;
    //
    glUseProgram(atom_program);
    glBindVertexArray(atom_vao);
    GLuint offLoc = glGetUniformLocation(atom_program, "offset");
    for(int x=0;x!=mult[0];++x){
        for(int y=0;y!=mult[1];++y){
            for(int z=0;z!=mult[2];++z){
                glUniform3fv(offLoc, 1, (-center+x*xvec+y*yvec+z*zvec).data());
                glDrawArraysInstanced(GL_TRIANGLES,0,atom_model_npoly,atom_buffer.size());
            }
        }
    }
}

void GLWidget::resizeGL(int w, int h)
{
    h==0?h=1:0;
    glViewport(0,0,w,h);
    float aspect = float(w)/h;
    pMat = guiMatMkOrtho(-10*aspect, 10*aspect, -10, 10, 0, 1000);
    pMatChanged = true;;
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

void GLWidget::setStep(const Step* step)
{
    curStep = step;
    //atoms
    const std::vector<Atom>& atoms = step->getAtoms();
    atom_buffer.reserve(step->getNat());
    atom_buffer.clear();
    for(const Atom& at:atoms){
        PseEntry &pse = (*step->pse)[at.name];
        atom_buffer.push_back({{at.coord[0],at.coord[1],at.coord[2],pse.covr,
                          pse.col[0],pse.col[1],pse.col[2],pse.col[3]}});
        atoms_changed = true;
    }
    update();
}

void GLWidget::setCamera(int i)
{
    switch(i){
    case -2: // +x
        rMat = {{ 0, 1, 0, 0,
                      0, 0, 1, 0,
                      1, 0, 0, 0,
                      0, 0, 0, 1}};
        break;
    case -5: // -x
        rMat = {{ 0,-1, 0, 0,
                      0, 0, 1, 0,
                     -1, 0, 0, 0,
                      0, 0, 0, 1}};
        break;
    case -3: // +y
        rMat = {{-1, 0, 0, 0,
                      0, 0, 1, 0,
                      0, 1, 0, 0,
                      0, 0, 0, 1}};
        break;
    case -6: // -y
        rMat = {{ 1, 0, 0, 0,
                      0, 0, 1, 0,
                      0,-1, 0, 0,
                      0, 0, 0, 1}};
        break;
    case -4: // +z
        rMat = {{ 1, 0, 0, 0,
                      0, 1, 0, 0,
                      0, 0, 1, 0,
                      0, 0, 0, 1}};
        break;
    case -7: // -z
        rMat = {{-1, 0, 0, 0,
                      0, 1, 0, 0,
                      0, 0,-1, 0,
                      0, 0, 0, 1}};
        break;
    }
    rMatChanged = true;
    update();
}

void GLWidget::keyPressEvent(QKeyEvent *e)
{
    switch(e->key()){
        case Qt::Key_Down:
            guiMatRot(rMat, 10, 1, 0, 0);
            rMatChanged = true;
            break;
        case Qt::Key_Up:
            guiMatRot(rMat, -10, 1, 0, 0);
            rMatChanged = true;
            break;
        case Qt::Key_Left:
            guiMatRot(rMat, -10, 0, 1, 0);
            rMatChanged = true;
            break;
        case Qt::Key_Right:
            guiMatRot(rMat, 10, 0, 1, 0);
            rMatChanged = true;
            break;
        default:
            return;
    }
    e->accept();
    update();
}

void GLWidget::wheelEvent(QWheelEvent *e)
{
    guiMatScale(vMat, e->angleDelta().y()>0?1.1:0.9);
    vMatChanged = true;
    e->accept();
    update();
}

void mouseEvent(){

}

void GLWidget::mousePressEvent(QMouseEvent *e)
{
    setFocus();
    if (!(e->buttons()&7)) return;
    mousePos = e->pos();
    switch(mouseMode){
        case MouseMode::Camera:
            if(e->button() == Qt::MouseButton::RightButton){
                rMat = {{1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1}};
                rMatChanged = true;
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
                guiMatRot(rMat, delta.x(), 0, 1, 0);
                guiMatRot(rMat, delta.y(), 1, 0, 0);
                rMatChanged = true;
                update();
            }else if(e->buttons() & Qt::MouseButton::MiddleButton){
                guiMatTranslate(vMat, delta.x()/10., -delta.y()/10., 0);
                vMatChanged = true;
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
