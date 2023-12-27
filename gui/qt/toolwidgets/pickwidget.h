#ifndef PICKWIDGET_H
#define PICKWIDGET_H

#include <QWidget>

#include "vipster/molecule.h"

namespace Ui {
class PickWidget;
}

class PickWidget : public QWidget
{
    Q_OBJECT

public:
    explicit PickWidget(QWidget *parent = nullptr);
    ~PickWidget() override;
    void updateSelection(const Vipster::Step::selection &);

private:
    Ui::PickWidget *ui;
};

#endif // PICKWIDGET_H
