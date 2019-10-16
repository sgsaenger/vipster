#ifndef PYTHONWIDGET_PY_H
#define PYTHONWIDGET_PY_H

#include <QWidget>
#include "../basewidget.h"

namespace Ui {
class PythonWidget;
}

class PythonWidget : public BaseWidget
{
    Q_OBJECT

public:
    explicit PythonWidget(QWidget *parent = nullptr);
    ~PythonWidget();

private:
    Ui::PythonWidget *ui;
};

#endif // PYTHONWIDGET_PY_H
