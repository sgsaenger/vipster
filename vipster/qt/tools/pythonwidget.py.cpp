#include "pythonwidget.py.h"
#include "ui_pythonwidget.py.h"

using namespace Vipster;

PythonWidget::PythonWidget(QWidget *parent) :
    BaseWidget(parent),
    ui(new Ui::PythonWidget)
{
    ui->setupUi(this);
}

PythonWidget::~PythonWidget()
{
    delete ui;
}
