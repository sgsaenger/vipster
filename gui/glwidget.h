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
#include "libvipster.h"

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
    //void mousePressEvent(QMouseEvent *e);
    //void mouseMoveEvent(QMouseEvent *e);
    //void mouseReleaseEvent(QMouseEvent *);
    void setStep(const Vipster::Step& step);
public slots:
    void setMode(int i,bool t);
private:
    QOpenGLVertexArrayObject vao;
    QOpenGLShaderProgram atom_shader,bond_shader,cell_shader;
    QOpenGLBuffer sphere_vbo,torus_vbo;     //model-geometries
    QOpenGLBuffer atom_vbo;                 //positions and properties
    std::array<QOpenGLBuffer,8> bond_vbo;   //gpu-side data
    QOpenGLBuffer cell_vbo; //gpu-side data
    QOpenGLBuffer cell_ibo;
    std::vector<std::array<float,8>> atom_buffer; //cpu-side data
    std::array<std::vector<std::array<float,24>>,8> bond_buffer; //cpu-side data
    std::array<std::array<float,3>,8> cell_buffer = {};
    QMatrix4x4 pMatrix,vMatrix,rMatrix;
    std::array<int,3> mult = {{1,1,1}}; //number of repetitions
    bool aVo,bVo,cVo; //keep track if vbos are outdated
    void drawAtoms(void);
    void drawBonds(void);
    void drawCell(void);
};

#endif // GLWIDGET_
