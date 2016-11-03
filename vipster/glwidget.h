#ifndef GLWIDGET_H
#define GLWIDGET_H

#include <QOpenGLFunctions_3_3_Core>
#include <QOpenGLWidget>
#include <QOpenGLShaderProgram>
#include <QOpenGLBuffer>
#include <QOpenGLVertexArrayObject>
#include <QApplication>
#include <iostream>
#include <vector>
#include <array>
#include <step.h>

class GLWidget: public QOpenGLWidget, protected QOpenGLFunctions_3_3_Core
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
    const Vipster::Step* curStep{NULL}; // Pointer to currently loaded Step
    // OGL-CPU/GPU buffers
    QOpenGLVertexArrayObject vao;
    QOpenGLShaderProgram atom_shader,bond_shader,cell_shader;
    QOpenGLBuffer sphere_vbo,torus_vbo;     //model-geometries
    QOpenGLBuffer atom_vbo;                 //positions and properties
    QOpenGLBuffer bond_vbo; //gpu-side data
    QOpenGLBuffer cell_vbo; //gpu-side data
    QOpenGLBuffer cell_ibo{QOpenGLBuffer::IndexBuffer}; //gpu-side data
    std::vector<std::array<float,8>> atom_buffer;  //cpu-side data
    std::vector<std::array<float,24>> bond_buffer; //cpu-side data
    std::array<std::array<float,3>,8> cell_buffer; //cpu-side data
    // Other data for rendering
    QMatrix4x4 pMatrix,vMatrix,rMatrix;
    std::array<int,3> mult{{1,1,1}}; //number of repetitions
    bool aVo{true},bVo{true},cVo{true}; //keep track if vbos are outdated
    float xshift{0.0}, yshift{0.0}, distance{1.0};
    // Input handling
    enum class MouseMode { Camera, Select, Modify };
    MouseMode mouseMode{MouseMode::Camera};
    QPoint mousePos;
    // private functions:
    void drawAtoms(void);
    void drawBonds(void);
    void drawCell(void);
};

#endif // GLWIDGET_
