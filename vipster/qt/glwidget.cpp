#include "glwidget.h"
#include <QKeyEvent>
#include "../common/atom_model.h"
#include "../common/bond_model.h"

using namespace Vipster;

GLWidget::GLWidget(QWidget *parent):
    QOpenGLWidget(parent) {}

GLWidget::~GLWidget()
{
    makeCurrent();

    doneCurrent();
}

void GLWidget::initializeGL()
{
    initializeOpenGLFunctions();
    glClearColor(1,1,1,1);
    glEnable(GL_MULTISAMPLE);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    initShaders("# version 330\n", ":/shaders");
    //
    initAtomVAO();
    initBondVAO();
    initCellVAO();
    initViewUBO();
    initViewMat();
}

void GLWidget::paintGL()
{
    updateViewUBO();
    //
    updateVBOs();
    //
    draw();
}

void GLWidget::resizeGL(int w, int h)
{
    resizeViewMat(w, h);
}

void GLWidget::setMode(int, bool t)
{
    // TODO: implement mouse-modes
    if(!t) {
        return;
    }
}

void GLWidget::setMult(int i)
{
    if(QObject::sender()->objectName() == "xMultBox"){ mult[0] = static_cast<uint8_t>(i); }
    else if(QObject::sender()->objectName() == "yMultBox"){ mult[1] = static_cast<uint8_t>(i); }
    else if(QObject::sender()->objectName() == "zMultBox"){ mult[2] = static_cast<uint8_t>(i); }
    update();
}

void GLWidget::setStep(const StepProper* step)
{
    updateBuffers(step, settings.showBonds.val);
    update();
}

void GLWidget::setCamera(int i)
{
    alignViewMat(static_cast<alignDir>((-i)-2));
    update();
}

void GLWidget::keyPressEvent(QKeyEvent *e)
{
    switch(e->key()){
    case Qt::Key_Down:
        rotateViewMat(0, 10, 0);
        break;
    case Qt::Key_Up:
        rotateViewMat(0, -10, 0);
        break;
    case Qt::Key_Left:
        rotateViewMat(-10, 0, 0);
        break;
    case Qt::Key_Right:
        rotateViewMat(10, 0, 0);
        break;
    default:
        return;
    }
    e->accept();
    update();
}

void GLWidget::wheelEvent(QWheelEvent *e)
{
    zoomViewMat(e->angleDelta().y());
    e->accept();
    update();
}

void GLWidget::mousePressEvent(QMouseEvent *e)
{
    setFocus();
    if (!(e->buttons()&7)) {
        return;
    }
    mousePos = e->pos();
    switch(mouseMode){
        case MouseMode::Camera:
            if(e->button() == Qt::MouseButton::RightButton){
                alignViewMat(alignDir::z);
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
    if (!(e->buttons()&7)) {
        return;
    }
    QPoint delta = e->pos() - mousePos;
    switch(mouseMode){
        case MouseMode::Camera:
            if((e->buttons() & Qt::MouseButton::LeftButton) != 0u){
                rotateViewMat(delta.x(), delta.y(), 0);
                update();
            }else if((e->buttons() & Qt::MouseButton::MiddleButton) != 0u){
                translateViewMat(delta.x(), -delta.y(), 0);
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
    if (!(e->buttons()&7)) {
        return;
    }
    switch(mouseMode){
        case MouseMode::Camera:
            break;
        case MouseMode::Select:
            break;
        case MouseMode::Modify:
            break;
    }
}
