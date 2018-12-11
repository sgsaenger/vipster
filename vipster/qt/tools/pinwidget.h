#ifndef PINWIDGET_H
#define PINWIDGET_H

#include <QWidget>
#include "../mainwindow.h"

namespace Ui {
class PinWidget;
}

class PinWidget : public BaseWidget
{
    Q_OBJECT

public:
    explicit PinWidget(QWidget *parent = nullptr);
    ~PinWidget() override;
    void updateWidget(Vipster::guiChange_t) override;

private slots:
    void on_showStep_toggled(bool checked);

    void on_showCell_toggled(bool checked);

    void on_delStep_clicked();

    void on_addStep_clicked();

    void on_stepList_currentRowChanged(int currentRow);

    void on_insertStep_clicked();

private:
    Ui::PinWidget *ui;
    struct PinnedStep{
        bool display, cell;
        Vipster::GUI::StepData gpu_data;
    };
    std::vector<Vipster::Step*> stepList;
    std::map<Vipster::Step*, PinnedStep> stepMap;
    Vipster::Step* mainStep{nullptr};
    Vipster::Step* activeStep{nullptr};
};

#endif // PINWIDGET_H
