#ifndef BONDWIDGET_H
#define BONDWIDGET_H

#include <QWidget>
#include "../mainwindow.h"
#include "bondmodel.h"

namespace Ui {
class BondWidget;
}

class BondWidget : public BaseWidget
{
    Q_OBJECT

public:
    explicit BondWidget(QWidget *parent = nullptr);
    ~BondWidget() override;
    void updateWidget(Vipster::GUI::change_t change) override;

private:
    Ui::BondWidget *ui;
    BondModel bondModel{this};
};

#endif // BONDWIDGET_H
