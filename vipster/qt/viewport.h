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
    explicit ViewPort(QWidget *parent = nullptr, bool active=false);
    explicit ViewPort(const ViewPort &vp);
    ~ViewPort() override;
    void triggerUpdate(Vipster::GUI::change_t change) override;
    void updateWidget(Vipster::GUI::change_t change) override;
    void registerMol(const std::string& name);
    void makeActive(bool active);

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
    void on_closeButton_clicked();
    void on_vSplitButton_clicked();
    void on_hSplitButton_clicked();
    void on_checkActive_toggled(bool checked);

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
