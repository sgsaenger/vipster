#ifndef VIEWPORT_H
#define VIEWPORT_H

#include <QAbstractButton>
#include <QTimer>
#include <QKeyEvent>
#include <QFrame>

#include "vipster/molecule.h"
#include "guiglobals.h"
#include "guidata.h"
#include "seldata.h"

namespace Ui {
class ViewPort;
}

class GLWidget;
class MainWindow;

class ViewPort : public QFrame
{
    Q_OBJECT

    friend class MainWindow;
    friend class GLWidget;
public:
    explicit ViewPort(MainWindow *parent);
    explicit ViewPort(const ViewPort &vp);
    ~ViewPort() override;

    void addExtraData(const std::shared_ptr<Vipster::GUI::Data> &dat, bool global);
    void delExtraData(const std::shared_ptr<Vipster::GUI::Data> &dat, bool global);
    bool hasExtraData(const std::shared_ptr<Vipster::GUI::Data> &dat, bool global);
    void cleanExtraData();

    // Extra visualizations attached to this viewport
    struct {
        std::vector<std::weak_ptr<Vipster::GUI::Data>> extras{};
    } vpdata{};

    // Visualization state of this molecule in this viewport
    struct MolState{
        size_t curStep{0};                  // current step
        Vipster::GUI::PBCVec mult{1,1,1};   // crystal multiplier
    };
    std::map<const Vipster::Molecule*, MolState> moldata{};

    // State of this step in this viewport
    struct StepState{
        std::unique_ptr<Vipster::Step::selection> sel{nullptr};     // selection
        std::vector<std::weak_ptr<Vipster::GUI::Data>> extras{};    // per step extra visualizations
    };
    std::map<const Vipster::Step*, StepState> stepdata{};

public slots:
    void setMol(int i);
    void setStep(int i, bool setMol=false);
    void setMult(int i);
    void setMultEnabled(bool b);
    void setBondMode(bool b);
    void stepBut(QAbstractButton *but);
    void setMouseMode(int i);
    void cameraBut(QAbstractButton *but);

private slots:
    void on_mouseMode_currentIndexChanged(int index);
    void updateMoleculeList(const std::list<Vipster::Molecule> &molecules);

public:
    void updateState();
    bool isActive();

    Ui::ViewPort *ui;
    MainWindow *master;
    GLWidget *openGLWidget{nullptr};
    Vipster::Molecule* curMol{nullptr};
    Vipster::Step* curStep{nullptr};
    Vipster::Step::selection* curSel{nullptr};
    QTimer playTimer{};
};

#endif // VIEWPORT_H
