#ifndef PINWIDGET_H
#define PINWIDGET_H

#include <QWidget>
#include "stepdata.h"
#include "vipster/molecule.h"

namespace Ui {
class PinWidget;
}

class PinWidget : public QWidget
{
    Q_OBJECT

public:
    explicit PinWidget(QWidget *parent = nullptr);
    ~PinWidget() override;

private slots:
    void addPin();
    void delPin();
    void toggleVisible(bool checked);
    void toggleCell(bool checked);
    void setMult(int);
    void setOffset(double);
    void toggleRepeat(bool checked);
    void selectPin(int currentRow);
    void insertStep();

private:
    Ui::PinWidget *ui;
    struct PinnedStep: public Vipster::GUI::StepData{
        PinnedStep(const Vipster::Step *step, const std::string& name,
                   Vipster::GUI::PBCVec mult);
        std::string name;
        Vipster::GUI::PBCVec mult;
        Vipster::Vec offset{};
        bool repeat{true}, cell{true};
        void draw(const Vipster::Vec &off, const Vipster::GUI::PBCVec &mult,
                  const Vipster::Mat &cv, bool drawCell, void *context) override;
    };
    // TODO: why shared_ptr instead of value?
    std::vector<std::shared_ptr<PinnedStep>> pinnedSteps;
    std::shared_ptr<PinnedStep> curPin{nullptr};
};

#endif // PINWIDGET_H
