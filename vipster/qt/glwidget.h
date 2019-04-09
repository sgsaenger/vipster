#ifndef GLWIDGET_H
#define GLWIDGET_H

#include <QOpenGLWidget>
#include "step.h"
#include "mainwindow.h"
#include "../common/guiwrapper.h"

class GLWidget: public QOpenGLWidget, public Vipster::GuiWrapper
{
    Q_OBJECT

public:
    explicit GLWidget(QWidget *parent = nullptr);
    ~GLWidget() override;
    void initializeGL(void) override;
    void paintGL(void) override;
    void resizeGL(int w, int h) override;
    void keyPressEvent(QKeyEvent *e) override;
    void wheelEvent(QWheelEvent *e) override;
    void mousePressEvent(QMouseEvent *e) override;
    void mouseMoveEvent(QMouseEvent *e) override;
    void mouseReleaseEvent(QMouseEvent *e) override;
    void triggerUpdate(Vipster::guiChange_t change);
    void updateWidget(Vipster::guiChange_t change);
public slots:
    void setMode(int i,bool t);
    void setMult(int i);
    void setCamera(int i);
private:
    bool updateTriggered{false};
    MainWindow* master;
    // Input handling
    enum class MouseMode { Camera=-2, Select=-3, Modify=-4 };
    MouseMode mouseMode{MouseMode::Camera};
    QPoint mousePos, rectPos;
    Vipster::Vec shift;
    std::map<size_t, std::vector<Vipster::SizeVec>> pickAtoms();
    void rotAtoms(QPoint delta);
    void shiftAtomsXY(QPoint delta);
    void shiftAtomsZ(QPoint delta);
};

#endif // GLWIDGET_
