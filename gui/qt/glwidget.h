#ifndef GLWIDGET_H
#define GLWIDGET_H

#include <QOpenGLWidget>
#include "vipster/step.h"
#include "viewport.h"
#include "guiwrapper.h"

class GLWidget: public QOpenGLWidget, public Vipster::GuiWrapper
{
    Q_OBJECT

public:
    enum class MouseMode { Camera, Select, Modify, Bond };
    explicit GLWidget(QWidget *parent, const Vipster::Settings &settings);
    ~GLWidget() override;
    void initializeGL(void) override;
    void paintGL(void) override;
    void resizeGL(int w, int h) override;
    void keyPressEvent(QKeyEvent *e) override;
    void wheelEvent(QWheelEvent *e) override;
    void mousePressEvent(QMouseEvent *e) override;
    void mouseMoveEvent(QMouseEvent *e) override;
    void mouseReleaseEvent(QMouseEvent *e) override;
    void focusInEvent(QFocusEvent *e) override;
    void triggerUpdate(Vipster::GUI::change_t change);
    void updateWidget(Vipster::GUI::change_t change);
    void setMult(Vipster::GUI::PBCVec mult);
    void setMouseMode(MouseMode i);
    void setCamera(int i);
private:
    bool updateTriggered{false};
    // Input handling
    MouseMode mouseMode{MouseMode::Camera};
    QPoint mousePos, rectPos;
    Vipster::Vec shift;
    Vipster::SelectionIndices pickAtoms(QPoint from, QPoint to);
    void rotAtoms(QPoint delta);
    void shiftAtomsXY(QPoint delta);
    void shiftAtomsZ(QPoint delta);
};

#endif // GLWIDGET_
