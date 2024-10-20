#include "pythonwidget.py.h"
#include "ui_pythonwidget.py.h"

PythonWidget::PythonWidget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::PythonWidget)
{
    ui->setupUi(this);
}

PythonWidget::~PythonWidget()
{
    delete ui;
}
