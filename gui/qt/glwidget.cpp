#include "glwidget.h"
#include "mainwindow.h"
#include <QKeyEvent>
#include <QOpenGLFramebufferObject>
#include <QApplication>
#include <QPainter>

using namespace Vipster;

GLWidget::GLWidget(QWidget *parent, const Vipster::Settings& settings):
    QOpenGLWidget(parent),
    GuiWrapper{settings}
{
    if (QSurfaceFormat::defaultFormat().renderableType() == QSurfaceFormat::OpenGL){
        // set surface format required for the hacky selection algorithm
        setTextureFormat(GL_RGBA16F);
    }
    setFocusPolicy(Qt::WheelFocus);
    vpExtras = &static_cast<ViewPort*>(parent)->vpdata.extras;
}

GLWidget::~GLWidget()
{
    makeCurrent();

    doneCurrent();
}

void GLWidget::triggerUpdate(GUI::change_t change){
    updateTriggered = true;
    static_cast<ViewPort*>(parent())->triggerUpdate(change);
}

void GLWidget::updateWidget(GUI::change_t change)
{
    if((change & GUI::stepChanged) == GUI::stepChanged ){
        setMainStep(static_cast<ViewPort*>(parent())->curStep);
        setMainSel(static_cast<ViewPort*>(parent())->curSel);
        stepExtras = &static_cast<ViewPort*>(parent())->stepdata[curStep].extras;
    }else{
        if(change & GUI::Change::settings){
            selection.color = settings.selCol.val;
        }
        if(change & (GUI::Change::atoms | GUI::Change::cell | GUI::Change::fmt |
                     GUI::Change::settings)) {
            updateMainStep();
            updateMainSelection();
        }else if(change & GUI::Change::selection){
            updateMainSelection();
        }
    }
    update();
}

void GLWidget::initializeGL()
{
    makeCurrent();
    initGL();
}

void GLWidget::paintGL()
{
    QPainter painter{this};
    painter.beginNativePainting();
    draw(this);
    if(rectPos != mousePos){
        glDisable(GL_CULL_FACE);
        glDisable(GL_DEPTH_TEST);
        painter.endNativePainting();
        QPen pen{};
        pen.setWidth(2);
        pen.setStyle(Qt::SolidLine);
        pen.setColor(Qt::black);
        painter.setPen(pen);
        QBrush brush{};
        brush.setColor(QColor{180,180,180,40});
        brush.setStyle(Qt::SolidPattern);
        painter.setBrush(brush);
        painter.drawRect(QRect{mousePos, rectPos});
    }else{
        painter.endNativePainting();
    }
}

void GLWidget::resizeGL(int w, int h)
{
    resizeViewMat(w, h);
}

void GLWidget::setMouseMode(MouseMode mode)
{
    mouseMode = mode;
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
    case MouseMode::Bond:
        setCursor(Qt::CrossCursor);
        break;
    }
}

void GLWidget::setMult(GUI::PBCVec m)
{
    mult = m;
    update();
}

void GLWidget::setCamera(alignDir dir)
{
    alignViewMat(dir);
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
    case Qt::Key_R:
        static_cast<ViewPort*>(parent())->setMouseMode(0);
        break;
    case Qt::Key_S:
        static_cast<ViewPort*>(parent())->setMouseMode(1);
        break;
    case Qt::Key_M:
        static_cast<ViewPort*>(parent())->setMouseMode(2);
        break;
    case Qt::Key_B:
        static_cast<ViewPort*>(parent())->setMouseMode(3);
        break;
    default:
        return QWidget::keyPressEvent(e);
    }
    e->accept();
    update();
}

SelectionIndices GLWidget::pickAtoms(QPoint from, QPoint to)
{
    // TODO: limited to 2^16 atoms per unit cell & 2^16 unit cells
    // return indices of atoms enclosed by rectangle
    // defined by `from` and `to`
    std::set<SelectionPair> idx;
    makeCurrent();
    drawSel(this);
    QOpenGLFramebufferObjectFormat format;
    format.setSamples(0);
    format.setInternalTextureFormat(GL_RGBA16F);
    // TODO: benchmark smaller framebuffer
//    QOpenGLFramebufferObject fbo{QSize{w, h}, format};
    QOpenGLFramebufferObject fbo{size(), format};
    glBindFramebuffer(GL_READ_FRAMEBUFFER, defaultFramebufferObject());
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, fbo.handle());
    glBlitFramebuffer(0, 0, width(), height(),
                      0, 0, fbo.width(), fbo.height(),
                      GL_COLOR_BUFFER_BIT, GL_NEAREST);
    fbo.bind();
    auto x = std::min(from.x(), to.x());
    auto y = std::min(height() - 1 - from.y(),
                      height() - 1 - to.y());
    auto w = std::max(1, std::abs(from.x() - to.x()));
    auto h = std::max(1, std::abs(from.y() - to.y()));
    std::vector<GLuint64> data(static_cast<size_t>(w*h));
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glPixelStorei(GL_PACK_ALIGNMENT, 2);
    glReadPixels(x, y, w, h, GL_RGBA, GL_UNSIGNED_SHORT, data.data());
    fbo.release();
    std::set<GLuint64> set{data.begin(), data.end()};
    //
    for(const auto& p: set){
        auto i = static_cast<GLuint>(p&0xFFFFFFFF);
        auto box = static_cast<GLuint>((p&0xFFFFFFFF00000000)>>32);
        if(box){
            //untangle pbc-id
            box -= 1;
            auto z = box / (mult[0]*mult[1]);
            auto yrem = box % (mult[0]*mult[1]);
            auto y = yrem / mult[0];
            auto x = yrem % mult[0];
            // untangle atom-id
            idx.insert({i, {x,y,z}});
        }
    }
    return {idx.begin(), idx.end()};
}

void GLWidget::rotAtoms(QPoint delta)
{
    if(delta.isNull()){
        return;
    }
    double angle = delta.manhattanLength();
    auto axes = getAxes();
    Vec axis = delta.y() * axes[0] + -delta.x() * axes[1];
    if(curSel->getNat()){
        curSel->asFmt(AtomFmt::Bohr).modRotate(angle, axis, shift);
    }else{
        curStep->asFmt(AtomFmt::Bohr).modRotate(angle, axis, shift);
    }
    triggerUpdate(GUI::Change::atoms);
}

void GLWidget::shiftAtomsXY(QPoint delta)
{
    auto axes = Mat_trans(Mat_inv(getAxes()));
    Vec axis = delta.x() * axes[0] + delta.y() * axes[1];
    if(curSel->getNat()){
        curSel->asFmt(AtomFmt::Bohr).modShift(axis, 0.01);
    }else{
        curStep->asFmt(AtomFmt::Bohr).modShift(axis, 0.01);
    }
    triggerUpdate(GUI::Change::atoms);
}

void GLWidget::shiftAtomsZ(QPoint delta)
{
    auto axes = Mat_trans(Mat_inv(getAxes()));
    double fac = 0.01 * (delta.x() + delta.y());
    if(curSel->getNat()){
        curSel->asFmt(AtomFmt::Bohr).modShift(axes[2], fac);
    }else{
        curStep->asFmt(AtomFmt::Bohr).modShift(axes[2], fac);
    }
    triggerUpdate(GUI::Change::atoms);
}

void GLWidget::wheelEvent(QWheelEvent *e)
{
    zoomViewMat(e->angleDelta().y()>0?1.1f:0.9f);
    e->accept();
    update();
}

void GLWidget::mousePressEvent(QMouseEvent *e)
{
    e->accept();
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
            *curSel = curStep->select(SelectionFilter{});
            triggerUpdate(GUI::Change::selection);
        }
        break;
    case MouseMode::Modify:
        if((e->button() == Qt::MouseButton::LeftButton) &&
           (e->modifiers() & (Qt::Modifier::CTRL|Qt::Modifier::SHIFT)) == 0u){
            auto idx = pickAtoms(mousePos, mousePos);
            if(idx.empty()){
                if(curSel->getNat()){
                    shift = curSel->getCom(AtomFmt::Bohr);
                }else{
                    shift = curStep->getCom(AtomFmt::Bohr);
                }
            }else{
                SelectionFilter filter{};
                filter.mode = SelectionFilter::Mode::Index;
                filter.indices.push_back(*idx.begin());
                shift = curStep->asFmt(AtomFmt::Bohr).select(filter)[0].coord;
            }
        }
        break;
    case MouseMode::Bond:
        // nop
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
            if(e->modifiers() & Qt::Modifier::SHIFT){
                translateViewMat(delta.x(), -delta.y(), 0);
            }else{
                rotateViewMat(delta.x(), delta.y(), 0);
            }
            update();
        }else if((e->buttons() & Qt::MouseButton::MiddleButton) != 0u){
            translateViewMat(delta.x(), -delta.y(), 0);
            update();
        }
        rectPos = mousePos = e->pos();
        break;
    case MouseMode::Select:
        if((e->buttons() & Qt::MouseButton::RightButton) == 0u &&
                delta.manhattanLength() > 5){
            rectPos = e->pos();
            update();
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
        rectPos = mousePos = e->pos();
        break;
    case MouseMode::Bond:
        // nop
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
            // get selected atoms
            auto pick = pickAtoms(mousePos, rectPos);
            // join with existing selection if requested
            bool add = (e->button() & Qt::MouseButton::MiddleButton) ||
                       (e->modifiers() & Qt::ControlModifier);
            if(add){
                auto idx = curSel->getAtoms().indices;
                for(const auto& p: pick){
                    auto pos = std::find_if(idx.begin(), idx.end(), [&p](const auto& i){
                        return i == p;
                    });
                    if(pos == idx.end()){
                        idx.push_back(p);
                    }
                }
                std::swap(pick, idx);
            }
            // create new filter
            SelectionFilter filter{};
            filter.mode = SelectionFilter::Mode::Index;
            filter.indices = std::move(pick);
            // create new selection from filter
            *curSel = curStep->select(filter);
            rectPos = mousePos;
            triggerUpdate(GUI::Change::selection);
        }
        break;
    case MouseMode::Modify:
        //TODO: undo
        break;
    case MouseMode::Bond:
        // if we have manual bonds and picked exactly two atoms, toggle bond
        auto idx1 = pickAtoms(rectPos, rectPos);
        if(idx1.size() != 1) break;
        auto idx2 = pickAtoms(e->pos(), e->pos());
        if(idx2.size() != 1) break;
        auto at1 = *idx1.begin();
        auto at2 = *idx2.begin();
        DiffVec off_l = {static_cast<DiffVec::value_type>(at1.second[0]-at2.second[0]),
                         static_cast<DiffVec::value_type>(at1.second[1]-at2.second[1]),
                         static_cast<DiffVec::value_type>(at1.second[2]-at2.second[2]),
                        };
        DiffVec off_r = {-off_l[0], -off_l[1], -off_l[2]};
        // ignore bonds between one atom with itself
        if((at1.first == at2.first) && (off_l == DiffVec{0,0,0}))
            break;
        const auto& bonds = curStep->getBonds();
        // if bond is already present, delete it and return
        for(size_t i=0; i<bonds.size(); ++i){
            const auto& b = bonds[i];
            if(((b.at1 == at1.first) && (b.at2 == at2.first) && (b.diff == off_r)) ||
               ((b.at1 == at2.first) && (b.at2 == at1.first) && (b.diff == off_l))){
                curStep->delBond(i);
                triggerUpdate(GUI::Change::atoms);
                return;
            }
        }
        // create new bond
        curStep->addBond(at1.first, at2.first, off_r);
        triggerUpdate(GUI::Change::atoms);
        break;
    }
}

void GLWidget::focusInEvent(QFocusEvent *)
{
    // make sure our viewport is the active one
    auto vp = static_cast<ViewPort*>(parentWidget());
    if(!vp->active){
        vp->master->changeViewports(vp, MainWindow::VP_ACTIVE);
    }
}
