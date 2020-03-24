#ifndef PINWIDGET_H
#define PINWIDGET_H

#include <QWidget>
#include "../basewidget.h"
#include "stepdata.h"
#include "molecule.h"

namespace Ui {
class PinWidget;
}

class PinWidget : public BaseWidget
{
    Q_OBJECT

public:
    explicit PinWidget(QWidget *parent = nullptr);
    ~PinWidget() override;
    void updateWidget(Vipster::GUI::change_t) override;

private slots:
    void on_showCell_toggled(bool checked);
    void on_repeatStep_toggled(bool checked);
    void on_delStep_clicked();
    void on_addStep_clicked();
    void on_stepList_currentRowChanged(int currentRow);
    void on_insertStep_clicked();
    void setMult(int);
    void setOffset(double);
    void on_showStep_toggled(bool checked);

    void on_helpButton_clicked();

private:
    Ui::PinWidget *ui;
    struct PinnedStep: public Vipster::GUI::StepData{
        PinnedStep(const Vipster::GUI::GlobalData &glob, Vipster::Step *step,
                   const std::string& name, Vipster::GUI::PBCVec mult);
        std::string name;
        Vipster::GUI::PBCVec mult;
        Vipster::Vec offset{};
        bool repeat{true}, cell{true};
        void draw(const Vipster::Vec &off, const Vipster::GUI::PBCVec &mult,
                  const Vipster::Mat &cv, bool drawCell, void *context) override;
    };
    std::vector<std::shared_ptr<PinnedStep>> pinnedSteps;
    std::shared_ptr<PinnedStep> curPin{nullptr};
};

#endif // PINWIDGET_H
