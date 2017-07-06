#ifndef GLWIDGET_H
#define GLWIDGET_H

#include <QOpenGLWidget>
#include <QOpenGLShaderProgram>
#include <QOpenGLBuffer>
#include <QOpenGLVertexArrayObject>
#include <QApplication>
#include <iostream>
#include <vector>
#include <array>
#include <step.h>
#include <guiwrapper.h>

class GLWidget: public QOpenGLWidget
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
    void setStep(const Vipster::Step* step);
public slots:
    void setMode(int i,bool t);
    void setMult(int i);
    void setCamera(int i);
private:
    const Vipster::Step* curStep{nullptr}; // Pointer to currently loaded Step
    GuiWrapper gui;
    std::string readShader(QString filePath);
    // Other data for rendering
    std::array<int,3> mult{{1,1,1}}; //number of repetitions
    float xshift{0.0}, yshift{0.0}, distance{1.0};
    // Input handling
    enum class MouseMode { Camera, Select, Modify };
    MouseMode mouseMode{MouseMode::Camera};
    QPoint mousePos;
};

#endif // GLWIDGET_
