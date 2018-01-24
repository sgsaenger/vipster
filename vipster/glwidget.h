#ifndef GLWIDGET_H
#define GLWIDGET_H

#include <QOpenGLWidget>
#include "step.h"
#include "guiwrapper.h"

class GLWidget: public QOpenGLWidget, private Vipster::GuiWrapper
{
    Q_OBJECT

public:
    explicit GLWidget(QWidget *parent = 0);
    ~GLWidget();
    void initializeGL(void);
    void paintGL(void);
    void resizeGL(int w, int h);
    void keyPressEvent(QKeyEvent *e);
    void wheelEvent(QWheelEvent *e);
    void mousePressEvent(QMouseEvent *e);
    void mouseMoveEvent(QMouseEvent *e);
    void mouseReleaseEvent(QMouseEvent *e);
    void setStep(const Vipster::StepProper* step);
public slots:
    void setMode(int i,bool t);
    void setMult(int i);
    void setCamera(int i);
private:
    // Input handling
    enum class MouseMode { Camera, Select, Modify };
    MouseMode mouseMode{MouseMode::Camera};
    QPoint mousePos;
};

#endif // GLWIDGET_
