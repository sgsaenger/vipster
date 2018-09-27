#include "glwidget.h"
#include <QKeyEvent>
#include <QOpenGLFramebufferObject>
#include <QApplication>

using namespace Vipster;

GLWidget::GLWidget(QWidget *parent):
    QOpenGLWidget(parent)
{
    for(auto *w: qApp->topLevelWidgets()){
        if(auto *t = qobject_cast<MainWindow*>(w)){
            master = t;
            return;
        }
    }
    throw Error("Could not determine MainWindow-instance.");
}

GLWidget::~GLWidget()
{
    makeCurrent();

    doneCurrent();
}

void GLWidget::triggerUpdate(uint8_t change){
    updateTriggered = true;
    master->updateWidgets(change);
}

void GLWidget::updateWidget(uint8_t change)
{
    if((change & guiStepChanged) == guiStepChanged ){
        setMainStep(master->curStep, settings.showBonds.val);
        setMainSel(master->curSel);
    }else{
        if(change & (GuiChange::atoms | GuiChange::cell | GuiChange::settings)) {
            updateMainStep(settings.showBonds.val);
            updateMainSelection();
        }else if(change & GuiChange::selection){
            updateMainSelection();
        }
    }
    update();
}

void GLWidget::initializeGL()
{
    makeCurrent();
    initGL("# version 330\n", ":/shaders");
}

void GLWidget::paintGL()
{
    draw();
}

void GLWidget::resizeGL(int w, int h)
{
    resizeViewMat(w, h);
}

void GLWidget::setMode(int mode, bool t)
{
    if(!t) {
        return;
    }
    mouseMode = static_cast<MouseMode>(mode);
    switch(mouseMode){
    case MouseMode::Camera:
        setCursor(Qt::ArrowCursor);
        break;
    case MouseMode::Modify:
        setCursor(Qt::OpenHandCursor);
        break;
    case MouseMode::Select:
        setCursor(Qt::CrossCursor);
        break;
    }
}

void GLWidget::setMult(int i)
{
    if(QObject::sender()->objectName() == "xMultBox"){ mult[0] = static_cast<uint8_t>(i); }
    else if(QObject::sender()->objectName() == "yMultBox"){ mult[1] = static_cast<uint8_t>(i); }
    else if(QObject::sender()->objectName() == "zMultBox"){ mult[2] = static_cast<uint8_t>(i); }
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

std::set<size_t> GLWidget::pickAtoms()
{
        std::set<size_t> idx;
        makeCurrent();
        drawSel();
        QOpenGLFramebufferObjectFormat format;
        format.setSamples(0);
        QOpenGLFramebufferObject fbo{size(), format};
        glBindFramebuffer(GL_READ_FRAMEBUFFER, defaultFramebufferObject());
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, fbo.handle());
        glBlitFramebuffer(0, 0, width(), height(),
                          0, 0, fbo.width(), fbo.height(),
                          GL_COLOR_BUFFER_BIT, GL_NEAREST);
        fbo.bind();
        auto x = std::min(mousePos.x(), rectPos.x());
        auto y = std::min(height() - 1 - mousePos.y(),
                          height() - 1 - rectPos.y());
        auto w = std::max(1, std::abs(mousePos.x() - rectPos.x()));
        auto h = std::max(1, std::abs(mousePos.y() - rectPos.y()));
        std::vector<GLubyte> data(4*static_cast<size_t>(w*h));
        glReadPixels(x, y, w, h, GL_RGBA, GL_UNSIGNED_BYTE, data.data());
        fbo.release();
        for(size_t i=0; i<data.size(); i+=4){
            if(data[i+3]){
                idx.insert(data[i+0] +
                           (static_cast<size_t>(data[i+1])<<8) +
                           (static_cast<size_t>(data[i+2])<<16)
                        );
            }
        }
        return idx;
}

void GLWidget::rotAtoms(QPoint delta)
{
    if(delta.isNull()){
        return;
    }
    float angle = delta.manhattanLength();
    auto axes = getAxes();
    Vec axis = delta.y() * axes[0] + delta.x() * axes[1];
    if(curSel->getNat()){
        curSel->modRotate(angle, axis, shift);
    }else{
        curStep->modRotate(angle, axis, shift);
    }
    triggerUpdate(GuiChange::atoms);
}

void GLWidget::shiftAtomsXY(QPoint delta)
{
    auto axes = Mat_trans(Mat_inv(getAxes()));
    Vec axis = delta.x() * axes[0] + delta.y() * axes[1];
    if(curSel->getNat()){
        curSel->modShift(axis, 0.01f);
    }else{
        curStep->modShift(axis, 0.01f);
    }
    triggerUpdate(GuiChange::atoms);
}

void GLWidget::shiftAtomsZ(QPoint delta)
{
    auto axes = Mat_trans(Mat_inv(getAxes()));
    float fac = 0.01f * (delta.x() + delta.y());
    if(curSel->getNat()){
        curSel->modShift(axes[2], fac);
    }else{
        curStep->modShift(axes[2], fac);
    }
    triggerUpdate(GuiChange::atoms);
}

void GLWidget::wheelEvent(QWheelEvent *e)
{
    zoomViewMat(e->angleDelta().y());
    e->accept();
    update();
}

void GLWidget::mousePressEvent(QMouseEvent *e)
{
    e->accept();
    setFocus();
    if (!(e->buttons()&7)) {
        return;
    }
    rectPos = mousePos = e->pos();
    switch(mouseMode){
    case MouseMode::Camera:
        if(e->button() == Qt::MouseButton::RightButton){
            alignViewMat(alignDir::z);
            update();
        }
        break;
    case MouseMode::Select:
        if(e->button() == Qt::MouseButton::RightButton){
            curSel->setFilter(SelectionFilter{});
            triggerUpdate(GuiChange::selection);
        }
        break;
    case MouseMode::Modify:
        if((e->button() == Qt::MouseButton::LeftButton) &&
           (e->modifiers() & (Qt::Modifier::CTRL|Qt::Modifier::SHIFT)) == 0u){
            auto idx = pickAtoms();
            if(idx.empty()){
                if(curSel->getNat()){
                    shift = curSel->getCom(curSel->getFmt());
                }else{
                    shift = curStep->getCom(curStep->getFmt());
                }
            }else{
                shift = (*curStep)[*idx.begin()].coord;
            }
        }
        break;
    }
}

void GLWidget::mouseMoveEvent(QMouseEvent *e)
{
    e->accept();
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
        mousePos = e->pos();
        break;
    case MouseMode::Select:
        if((e->buttons() & Qt::MouseButton::RightButton) == 0u &&
                delta.manhattanLength() > 5){
            rectPos = e->pos();
        }
        break;
    case MouseMode::Modify:
        switch(e->buttons()){
        case Qt::MouseButton::LeftButton:
            if(e->modifiers() & Qt::Modifier::CTRL){
                shiftAtomsXY(delta);
            }else if(e->modifiers() & Qt::Modifier::SHIFT){
                shiftAtomsZ(delta);
            }else{
                rotAtoms(delta);
            }
            break;
        case Qt::MouseButton::MiddleButton:
            shiftAtomsXY(delta);
            break;
        case Qt::MouseButton::RightButton:
            shiftAtomsZ(delta);
            break;
        default:
            break;
        }
        mousePos = e->pos();
        break;
    }
}

void GLWidget::mouseReleaseEvent(QMouseEvent *e)
{
    e->accept();
    switch(mouseMode){
    case MouseMode::Camera:
        break;
    case MouseMode::Select:
        if(e->button() != Qt::MouseButton::RightButton){
            bool add = (e->button() & Qt::MouseButton::MiddleButton) ||
                       (e->modifiers() & Qt::ControlModifier);
            SelectionFilter filter{};
            filter.mode = SelectionFilter::Mode::Index;
            if(add){
                const auto& origIndices = curSel->getIndices();
                filter.indices.insert(origIndices.begin(), origIndices.end());
            }
            auto idx = pickAtoms();
            filter.indices.insert(idx.begin(), idx.end());
            curSel->setFilter(filter);
            triggerUpdate(GuiChange::selection);
        }
        break;
    case MouseMode::Modify:
        //TODO: undo
        break;
    }
}
