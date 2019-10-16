#ifndef VIEWPORT_H
#define VIEWPORT_H

#include <QAbstractButton>
#include <QTimer>
#include <QKeyEvent>
#include <QFrame>

#include "molecule.h"
#include "toolwidgets.h"
#include "../common/guiglobals.h"

namespace Ui {
class ViewPort;
}

class GLWidget;
class MainWindow;

class ViewPort : public QFrame
{
    Q_OBJECT

    friend class MainWindow;
public:
    explicit ViewPort(MainWindow *parent, bool active=false);
    explicit ViewPort(const ViewPort &vp);
    ~ViewPort() override;
    void triggerUpdate(Vipster::GUI::change_t change);
    void updateWidget(Vipster::GUI::change_t change);
    void registerMol(const std::string& name);
    void makeActive(bool active);
    struct MolExtras{
        int curStep{-1};
        Vipster::GUI::PBCVec mult{1,1,1};
    };
    struct StepExtras{
        std::unique_ptr<Vipster::Step::selection> sel{nullptr};
        std::map<std::string, Vipster::Step::selection> def{};
    };
    std::map<Vipster::Molecule*, MolExtras> moldata;
    std::map<Vipster::Step*, StepExtras> stepdata;

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

public:
    Ui::ViewPort *ui;
    MainWindow *master;
    GLWidget *openGLWidget{nullptr};
    Vipster::Molecule* curMol{nullptr};
    Vipster::Step* curStep{nullptr};
    Vipster::Step::selection* curSel{nullptr};
    bool active{true};
    QTimer playTimer{};
};

#endif // VIEWPORT_H
