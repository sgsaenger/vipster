#ifndef MILLERWIDGET_H
#define MILLERWIDGET_H

#include <QWidget>
#include "../mainwindow.h"

namespace Ui {
class MillerWidget;
}

class MillerWidget : public BaseWidget
{
    Q_OBJECT

public:
    explicit MillerWidget(QWidget *parent = nullptr);
    ~MillerWidget() override;
    void updateWidget(uint8_t) override;

private slots:
    void updateIndex(int idx);
    void updateOffset(double off);
    void toggleView();

private:
    Ui::MillerWidget *ui;
    bool display{false};
    std::array<uint8_t, 3> plane;
    Vipster::Vec offset;
};

#endif // MILLERWIDGET_H
