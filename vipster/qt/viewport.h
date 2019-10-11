#ifndef VIEWPORT_H
#define VIEWPORT_H

#include <QAbstractButton>
#include <QTimer>
#include <QKeyEvent>

#include "basewidget.h"

#include "molecule.h"

namespace Ui {
class ViewPort;
}

class GLWidget;

class ViewPort : public BaseWidget
{
    Q_OBJECT

public:
    explicit ViewPort(QWidget *parent = nullptr);
    ~ViewPort();
    void updateWidget(Vipster::GUI::change_t) override;
    void registerMol(const std::string& name);

public slots:
    void setMol(int i);
    void setStep(int i, bool setMol=false);
    void setMult(int i);
    void setMultEnabled(bool b);
    void setBondMode(bool b);
    void stepBut(QAbstractButton *but);
    void setMouseMode(int i);
    void setCamera(int i);

private slots:
    void on_mouseMode_currentIndexChanged(int index);

private:
    friend class GLWidget;
    friend class MainWindow;
    Ui::ViewPort *ui;
    GLWidget *openGLWidget{nullptr};
    Vipster::Molecule* curMol{nullptr};
    Vipster::Step* curStep{nullptr};
    Vipster::Step::selection* curSel{nullptr};
    bool active{true};
    QTimer playTimer{};
};

#endif // VIEWPORT_H
