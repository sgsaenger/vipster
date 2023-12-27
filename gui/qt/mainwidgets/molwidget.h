#ifndef MOLWIDGET_H
#define MOLWIDGET_H

#include <QWidget>
#include <QItemSelection>
#include "../basewidget.h"

namespace Ui {
class MolWidget;
}

class MolWidget : public BaseWidget
{
    Q_OBJECT

public:
    explicit MolWidget(QWidget *parent = nullptr);
    ~MolWidget() override;

private:
    Ui::MolWidget *ui;
};

#endif // MOLWIDGET_H
