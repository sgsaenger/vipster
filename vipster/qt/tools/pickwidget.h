#ifndef PICKWIDGET_H
#define PICKWIDGET_H

#include <QWidget>
#include "../mainwindow.h"

namespace Ui {
class PickWidget;
}

class PickWidget : public BaseWidget
{
    Q_OBJECT

public:
    explicit PickWidget(QWidget *parent = nullptr);
    ~PickWidget() override;
    void updateWidget(Vipster::guiChange_t) override;

private:
    Ui::PickWidget *ui;
};

#endif // PICKWIDGET_H
